from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


if __name__ == "__main__":
    """
    Analyze the latest stats_*.csv produced by stats_snapping_distance_creation.py.
    """
    base_dir = Path(__file__).resolve().parent
    output_dir = base_dir.parent / "tests" / "output"
    candidates = sorted(output_dir.glob("stats_*.csv"))
    if not candidates:
        raise FileNotFoundError(
            f"No stats_*.csv found in {output_dir}. Run stats_snapping_distance_creation.py first."
        )

    csv_path = candidates[-1]
    df = pd.read_csv(csv_path, sep=";")
    print(f"using: {csv_path}")
    print(df.head())

    df["max"] = df.groupby(["key"])["diff"].transform("max")
    result = df.loc[df["diff"] == df["max"]].copy()
    result = result.sort_values(["distance", "key", "max"]).drop_duplicates(
        ["key", "max"]
    )
    print(result.head())

    fig1 = plt.figure()
    plt.scatter(result["max"], result["area"] / 1000)
    plt.xlabel("Max diff")
    plt.ylabel("Area / 1000")
    fig1.savefig(output_dir / "stats_scatter.png", dpi=150, bbox_inches="tight")

    ax = result.hist(
        column="distance",
        bins=50,
        grid=False,
        figsize=(10, 6),
        color="#86bf91",
    )
    fig2 = ax[0][0].get_figure()
    fig2.savefig(output_dir / "stats_distance_hist.png", dpi=150, bbox_inches="tight")

    print(f"plots written to: {output_dir}")
