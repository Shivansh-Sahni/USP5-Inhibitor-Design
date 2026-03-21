from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import (
    Descriptors,
    GraphDescriptors,
    Lipinski,
    MACCSkeys,
    rdFingerprintGenerator,
    rdMolDescriptors,
)
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

from usp5_workflow.data import annotate_labels, prepare_dataset, validate_and_canonicalize_smiles


class TanimotoKNN(BaseEstimator, RegressorMixin):
    def __init__(self, n_neighbors: int = 3, weighted: bool = True) -> None:
        self.n_neighbors = n_neighbors
        self.weighted = weighted

    def fit(self, X: list, y: np.ndarray) -> "TanimotoKNN":
        self.train_fps_ = list(X)
        self.y_ = np.asarray(y, dtype=float)
        return self

    def predict(self, X: list) -> np.ndarray:
        predictions = []
        for fp in X:
            sims = np.array([DataStructs.TanimotoSimilarity(fp, train_fp) for train_fp in self.train_fps_])
            idx = np.argsort(sims)[::-1][: self.n_neighbors]
            weights = sims[idx] + 1e-8
            values = self.y_[idx]
            if self.weighted:
                predictions.append(float(np.sum(weights * values) / np.sum(weights)))
            else:
                predictions.append(float(np.mean(values)))
        return np.array(predictions, dtype=float)


def _bitvector_to_array(bitvector) -> np.ndarray:
    arr = np.zeros((bitvector.GetNumBits(),), dtype=int)
    DataStructs.ConvertToNumpyArray(bitvector, arr)
    return arr


def load_variants() -> dict[str, pd.DataFrame]:
    raw = pd.read_csv("data/raw/First.csv")
    raw = validate_and_canonicalize_smiles(annotate_labels(raw))
    raw = raw[raw["is_valid_smiles"]].copy()
    raw["row_id"] = np.arange(len(raw))
    raw["pIC50"] = raw["pIC50"].astype(float)

    dedup = prepare_dataset(Path("data/raw/First.csv")).modeling.copy()
    dedup["row_id"] = np.arange(len(dedup))
    dedup["pIC50"] = dedup["target_pIC50"].astype(float)

    return {"raw_rows": raw, "dedup_rows": dedup}


def build_feature_sets(df: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], dict[str, list]]:
    base_rows = []
    graph_rows = []
    maccs_rows = []
    rdk_rows = []
    atom_pair_rows = []
    morgan_rows = []

    rdk_gen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=512)
    atom_pair_gen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=512)
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=256)

    fingerprints = {"morgan2_256": []}

    for _, row in df.iterrows():
        mol = row["mol"]
        row_id = row["row_id"]

        base_rows.append(
            {
                "row_id": row_id,
                "mw": Descriptors.MolWt(mol),
                "logp": Descriptors.MolLogP(mol),
                "tpsa": rdMolDescriptors.CalcTPSA(mol),
                "hbd": Lipinski.NumHDonors(mol),
                "hba": Lipinski.NumHAcceptors(mol),
                "rot": Lipinski.NumRotatableBonds(mol),
                "rings": Lipinski.RingCount(mol),
                "hac": Lipinski.HeavyAtomCount(mol),
                "fsp3": Lipinski.FractionCSP3(mol),
            }
        )
        graph_rows.append(
            {
                "row_id": row_id,
                "bertz": GraphDescriptors.BertzCT(mol),
                "balaban": GraphDescriptors.BalabanJ(mol),
                "chi0v": GraphDescriptors.Chi0v(mol),
                "chi1v": GraphDescriptors.Chi1v(mol),
                "chi2v": GraphDescriptors.Chi2v(mol),
                "kappa1": GraphDescriptors.Kappa1(mol),
                "kappa2": GraphDescriptors.Kappa2(mol),
                "kappa3": GraphDescriptors.Kappa3(mol),
            }
        )

        maccs = _bitvector_to_array(MACCSkeys.GenMACCSKeys(mol))
        maccs_rows.append({"row_id": row_id, **{f"maccs_{i}": int(v) for i, v in enumerate(maccs)}})

        rdk = _bitvector_to_array(rdk_gen.GetFingerprint(mol))
        rdk_rows.append({"row_id": row_id, **{f"rdk_{i}": int(v) for i, v in enumerate(rdk)}})

        atom_pair = _bitvector_to_array(atom_pair_gen.GetFingerprint(mol))
        atom_pair_rows.append({"row_id": row_id, **{f"ap_{i}": int(v) for i, v in enumerate(atom_pair)}})

        morgan_fp = morgan_gen.GetFingerprint(mol)
        fingerprints["morgan2_256"].append(morgan_fp)
        morgan = _bitvector_to_array(morgan_fp)
        morgan_rows.append({"row_id": row_id, **{f"m256_{i}": int(v) for i, v in enumerate(morgan)}})

    base_df = pd.DataFrame(base_rows)
    graph_df = pd.DataFrame(graph_rows)
    maccs_df = pd.DataFrame(maccs_rows)
    rdk_df = pd.DataFrame(rdk_rows)
    atom_pair_df = pd.DataFrame(atom_pair_rows)
    morgan_df = pd.DataFrame(morgan_rows)

    feature_sets = {
        "base": base_df,
        "graph": graph_df,
        "base_graph": base_df.merge(graph_df, on="row_id"),
        "maccs": maccs_df,
        "rdk512": rdk_df,
        "atom_pair512": atom_pair_df,
        "morgan2_256": morgan_df,
    }
    return feature_sets, fingerprints


def evaluate_feature_table(
    feature_df: pd.DataFrame,
    y: np.ndarray,
    model_name: str,
    model,
    variant_name: str,
    feature_name: str,
    evaluation: str,
) -> dict[str, object]:
    X = feature_df[[column for column in feature_df.columns if column != "row_id"]]

    if evaluation == "loocv":
        loo = LeaveOneOut()
        predictions = np.zeros(len(y), dtype=float)
        for train_idx, test_idx in loo.split(X):
            fitted = model.fit(X.iloc[train_idx], y[train_idx])
            predictions[test_idx] = float(np.ravel(fitted.predict(X.iloc[test_idx]))[0])
    else:
        predictions = np.ravel(model.fit(X, y).predict(X)).astype(float)

    return {
        "variant": variant_name,
        "feature_set": feature_name,
        "model": model_name,
        "evaluation": evaluation,
        "mae": mean_absolute_error(y, predictions),
        "rmse": math.sqrt(mean_squared_error(y, predictions)),
        "r2": r2_score(y, predictions),
    }


def main() -> None:
    variants = load_variants()
    output_path = Path("outputs/high_r2_search_results.csv")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, object]] = []

    standard_models = {
        "enet": Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("model", ElasticNet(alpha=0.01, l1_ratio=0.2, max_iter=30000, random_state=42)),
            ]
        ),
        "svr_rbf": Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("model", SVR(C=10.0, epsilon=0.1, kernel="rbf")),
            ]
        ),
        "rf": RandomForestRegressor(n_estimators=300, random_state=42, max_features="sqrt", min_samples_leaf=1),
        "extra": ExtraTreesRegressor(n_estimators=300, random_state=42, max_features="sqrt", min_samples_leaf=1),
    }
    trainfit_models = {
        "rf_trainfit": RandomForestRegressor(
            n_estimators=600, random_state=42, max_features="sqrt", min_samples_leaf=1
        ),
        "extra_trainfit": ExtraTreesRegressor(
            n_estimators=600, random_state=42, max_features="sqrt", min_samples_leaf=1
        ),
        "knn1_trainfit": Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                ("model", KNeighborsRegressor(n_neighbors=1)),
            ]
        ),
    }

    for variant_name, df in variants.items():
        feature_sets, fingerprints = build_feature_sets(df)
        y = df["pIC50"].to_numpy(dtype=float)

        for feature_name, feature_df in feature_sets.items():
            for model_name, model in standard_models.items():
                results.append(
                    evaluate_feature_table(
                        feature_df=feature_df,
                        y=y,
                        model_name=model_name,
                        model=model,
                        variant_name=variant_name,
                        feature_name=feature_name,
                        evaluation="loocv",
                    )
                )
            for model_name, model in trainfit_models.items():
                results.append(
                    evaluate_feature_table(
                        feature_df=feature_df,
                        y=y,
                        model_name=model_name,
                        model=model,
                        variant_name=variant_name,
                        feature_name=feature_name,
                        evaluation="train_fit",
                    )
                )

        for fp_name, fps in fingerprints.items():
            y_array = np.asarray(y, dtype=float)
            for model_name, model in {
                "tanimoto_knn1": TanimotoKNN(n_neighbors=1, weighted=False),
                "tanimoto_knn3w": TanimotoKNN(n_neighbors=3, weighted=True),
            }.items():
                loo = LeaveOneOut()
                predictions = np.zeros(len(y_array), dtype=float)
                for train_idx, test_idx in loo.split(np.arange(len(y_array))):
                    train_fps = [fps[i] for i in train_idx]
                    test_fps = [fps[i] for i in test_idx]
                    fitted = model.fit(train_fps, y_array[train_idx])
                    predictions[test_idx] = float(fitted.predict(test_fps)[0])
                results.append(
                    {
                        "variant": variant_name,
                        "feature_set": fp_name,
                        "model": model_name,
                        "evaluation": "loocv",
                        "mae": mean_absolute_error(y_array, predictions),
                        "rmse": math.sqrt(mean_squared_error(y_array, predictions)),
                        "r2": r2_score(y_array, predictions),
                    }
                )

            kernel = np.zeros((len(y_array), len(y_array)), dtype=float)
            for i in range(len(y_array)):
                for j in range(len(y_array)):
                    kernel[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])

            for alpha in [0.01, 0.1, 1.0]:
                model_name = f"tanimoto_krr_a{alpha:g}"
                loo = LeaveOneOut()
                predictions = np.zeros(len(y_array), dtype=float)
                for train_idx, test_idx in loo.split(kernel):
                    model = KernelRidge(alpha=alpha, kernel="precomputed")
                    model.fit(kernel[np.ix_(train_idx, train_idx)], y_array[train_idx])
                    predictions[test_idx] = float(model.predict(kernel[np.ix_(test_idx, train_idx)])[0])
                results.append(
                    {
                        "variant": variant_name,
                        "feature_set": fp_name,
                        "model": model_name,
                        "evaluation": "loocv",
                        "mae": mean_absolute_error(y_array, predictions),
                        "rmse": math.sqrt(mean_squared_error(y_array, predictions)),
                        "r2": r2_score(y_array, predictions),
                    }
                )

                model = KernelRidge(alpha=alpha, kernel="precomputed")
                model.fit(kernel, y_array)
                train_predictions = np.ravel(model.predict(kernel)).astype(float)
                results.append(
                    {
                        "variant": variant_name,
                        "feature_set": fp_name,
                        "model": f"{model_name}_trainfit",
                        "evaluation": "train_fit",
                        "mae": mean_absolute_error(y_array, train_predictions),
                        "rmse": math.sqrt(mean_squared_error(y_array, train_predictions)),
                        "r2": r2_score(y_array, train_predictions),
                    }
                )

    results_df = pd.DataFrame(results).sort_values(
        ["evaluation", "r2", "mae"], ascending=[True, False, True]
    )
    results_df.to_csv(output_path, index=False)

    honest_top = results_df[results_df["evaluation"] == "loocv"].head(10)
    optimistic_top = results_df[results_df["evaluation"] == "train_fit"].head(10)
    print("Top LOOCV results")
    print(honest_top.to_string(index=False))
    print("\nTop train-fit results")
    print(optimistic_top.to_string(index=False))


if __name__ == "__main__":
    main()
