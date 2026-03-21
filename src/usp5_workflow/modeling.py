from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.neighbors import KNeighborsRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

try:
    from xgboost import XGBRegressor
except Exception:  # pragma: no cover - optional dependency
    XGBRegressor = None


class MeanRegressor(BaseEstimator, RegressorMixin):
    def fit(self, X: pd.DataFrame, y: pd.Series) -> "MeanRegressor":
        self.mean_ = float(np.mean(y))
        return self

    def predict(self, X: pd.DataFrame) -> np.ndarray:
        return np.full(shape=(len(X),), fill_value=self.mean_)


@dataclass
class ModelResult:
    metrics: pd.DataFrame
    predictions: pd.DataFrame
    full_fit_predictions: pd.DataFrame


def _build_models(random_seed: int) -> dict[str, BaseEstimator]:
    models: dict[str, BaseEstimator] = {"mean": MeanRegressor()}

    for alpha in [0.01, 0.1, 0.5]:
        for l1_ratio in [0.2, 0.5, 0.8]:
            name = f"elastic_net_a{alpha:g}_l1{l1_ratio:g}"
            models[name] = Pipeline(
                steps=[
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    (
                        "model",
                        ElasticNet(
                            alpha=alpha,
                            l1_ratio=l1_ratio,
                            random_state=random_seed,
                            max_iter=20000,
                        ),
                    ),
                ]
            )

    for n_neighbors in [1, 3, 5]:
        for weights in ["uniform", "distance"]:
            name = f"knn_k{n_neighbors}_{weights}"
            models[name] = Pipeline(
                steps=[
                    ("imputer", SimpleImputer(strategy="median")),
                    ("scaler", StandardScaler()),
                    ("model", KNeighborsRegressor(n_neighbors=n_neighbors, weights=weights)),
                ]
            )

    for n_estimators in [200]:
        for min_samples_leaf in [1, 2]:
            for max_features in ["sqrt", 0.3]:
                for max_depth in [None, 4]:
                    depth_label = "none" if max_depth is None else str(max_depth)
                    feature_label = str(max_features).replace(".", "p")
                    name = (
                        f"random_forest_n{n_estimators}_leaf{min_samples_leaf}_"
                        f"feat{feature_label}_depth{depth_label}"
                    )
                    models[name] = RandomForestRegressor(
                        n_estimators=n_estimators,
                        random_state=random_seed,
                        min_samples_leaf=min_samples_leaf,
                        max_features=max_features,
                        max_depth=max_depth,
                    )

    if XGBRegressor is not None:
        for n_estimators in [100, 200]:
            for max_depth in [2, 3]:
                for learning_rate in [0.03, 0.1]:
                    name = f"xgboost_n{n_estimators}_d{max_depth}_lr{learning_rate:g}"
                    models[name] = XGBRegressor(
                        n_estimators=n_estimators,
                        max_depth=max_depth,
                        learning_rate=learning_rate,
                        subsample=0.9,
                        colsample_bytree=0.9,
                        reg_lambda=1.0,
                        objective="reg:squarederror",
                        random_state=random_seed,
                    )

    return models


def _metric_frame(model_name: str, feature_set: str, y_true: np.ndarray, y_pred: np.ndarray) -> pd.DataFrame:
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    try:
        r2 = r2_score(y_true, y_pred)
    except ValueError:
        r2 = np.nan

    return pd.DataFrame(
        [
            {
                "model": model_name,
                "feature_set": feature_set,
                "mae": mean_absolute_error(y_true, y_pred),
                "rmse": rmse,
                "r2": r2,
            }
        ]
    )


def run_loocv(
    modeling_df: pd.DataFrame,
    feature_sets: dict[str, pd.DataFrame],
    random_seed: int,
) -> ModelResult:
    metrics_frames = []
    prediction_frames = []
    full_fit_frames = []
    loo = LeaveOneOut()

    for feature_set_name, feature_df in feature_sets.items():
        merged = modeling_df.merge(feature_df, on="canonical_smiles", how="left")
        feature_columns = [column for column in feature_df.columns if column != "canonical_smiles"]
        X = merged[feature_columns]
        y = merged["target_pIC50"].to_numpy()

        for model_name, model in _build_models(random_seed).items():
            fold_predictions = np.zeros_like(y, dtype=float)
            for train_idx, test_idx in loo.split(X):
                X_train = X.iloc[train_idx]
                X_test = X.iloc[test_idx]
                y_train = y[train_idx]
                estimator = model
                estimator.fit(X_train, y_train)
                fold_predictions[test_idx] = estimator.predict(X_test)

            prediction_frame = merged[
                [
                    "representative_id",
                    "canonical_smiles",
                    "target_pIC50",
                    "row_count",
                    "target_origin",
                    "duplicate_group_has_conflict",
                    "n_measured_rows",
                    "n_assigned_rows",
                ]
            ].copy()
            prediction_frame["model"] = model_name
            prediction_frame["feature_set"] = feature_set_name
            prediction_frame["loocv_prediction"] = fold_predictions
            prediction_frame["residual"] = prediction_frame["target_pIC50"] - prediction_frame["loocv_prediction"]
            prediction_frames.append(prediction_frame)
            metrics_frames.append(_metric_frame(model_name, feature_set_name, y, fold_predictions))

            full_fit_model = model
            full_fit_model.fit(X, y)
            full_fit_frame = prediction_frame[
                ["representative_id", "canonical_smiles", "target_pIC50", "model", "feature_set"]
            ].copy()
            full_fit_frame["full_fit_prediction"] = full_fit_model.predict(X)
            full_fit_frames.append(full_fit_frame)

    metrics = pd.concat(metrics_frames, ignore_index=True).sort_values(["mae", "rmse", "model", "feature_set"])
    predictions = pd.concat(prediction_frames, ignore_index=True)
    full_fit_predictions = pd.concat(full_fit_frames, ignore_index=True)
    return ModelResult(metrics=metrics, predictions=predictions, full_fit_predictions=full_fit_predictions)


def build_ranked_compounds(
    modeling_df: pd.DataFrame,
    result: ModelResult,
) -> pd.DataFrame:
    best = result.metrics.iloc[0]
    ranked = result.full_fit_predictions[
        (result.full_fit_predictions["model"] == best["model"])
        & (result.full_fit_predictions["feature_set"] == best["feature_set"])
    ].copy()
    ranked = ranked.merge(
        modeling_df[
            [
                "canonical_smiles",
                "compound_ids",
                "target_origin",
                "row_count",
                "duplicate_group_has_conflict",
                "n_measured_rows",
                "n_assigned_rows",
            ]
        ],
        on="canonical_smiles",
        how="left",
    )
    ranked = ranked.sort_values("full_fit_prediction", ascending=False).reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    return ranked
