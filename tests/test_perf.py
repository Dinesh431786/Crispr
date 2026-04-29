import pandas as pd
import time
from crispr_app.analysis import find_off_targets_detailed

def test_off_target_performance():
    # Long random sequence
    import random
    bg = "".join(random.choice("ATCG") for _ in range(50000))
    # Add a near match
    guide = "ATGC" * 5
    bg = bg[:10000] + "ATGCATGCATGCATGCATGT" + "CGG" + bg[10023:]

    guides_df = pd.DataFrame({"gRNA": [guide]})

    start = time.time()
    results = find_off_targets_detailed(guides_df, bg, max_mismatches=3, pam="NGG")
    end = time.time()

    print(f"Time taken for 50kb scan: {end - start:.4f}s")
    assert not results.empty
    assert "gRNA" in results.columns

if __name__ == "__main__":
    test_off_target_performance()
