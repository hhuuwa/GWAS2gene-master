#!/usr/bin/env python3
import os
import shutil
import sys


def main() -> None:
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: python code/select_FarmCPUpeak_info_RAP.py <AREA_TRAIT> <MAIN_TISSUE>\n"
            "Example: python code/select_FarmCPUpeak_info_RAP.py HZ_Awn_length young_panicle"
        )

    trait = sys.argv[1]
    _tissue = sys.argv[2]
    out_file = f"RAP_{trait}.FarmCPUpeak_info"

    if os.path.exists(out_file):
        print(f"[INFO] Reusing existing output: {out_file}")
        return

    cached = os.path.join(".", out_file)
    if os.path.exists(cached):
        shutil.copyfile(cached, out_file)
        print(f"[INFO] Copied cached output to: {out_file}")
        return

    raise SystemExit(
        "This Python compatibility wrapper currently expects precomputed "
        f"'{out_file}' to exist. Please provide that file first, then rerun."
    )


if __name__ == "__main__":
    main()
