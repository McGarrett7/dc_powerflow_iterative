"""Download official MATPOWER IEEE benchmark cases into data/."""
from pathlib import Path
from urllib.request import urlopen

CASES = (14, 30, 57, 118, 300)
BASE_URL = "https://raw.githubusercontent.com/MATPOWER/matpower/master/data/case{case}.m"
DATA_DIR = Path(__file__).resolve().parent / "data"


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        target = DATA_DIR / f"case{case}.m"
        if target.exists():
            print(f"[SKIP] {target.name} đã tồn tại")
            continue
        url = BASE_URL.format(case=case)
        print(f"[GET] {url}")
        try:
            with urlopen(url, timeout=30) as response:
                content = response.read()
            if b"mpc.bus" not in content or b"mpc.branch" not in content:
                raise ValueError("Nội dung tải về không giống file MATPOWER case.")
            target.write_bytes(content)
            print(f"[OK] {target.name}: {len(content):,} bytes")
        except Exception as exc:
            print(f"[ERROR] Không tải được {target.name}: {exc}")
            print("        Có thể tải thủ công từ thư mục data/ của repository MATPOWER.")
    print("\nHoàn tất. Kiểm tra các file trong thư mục data/ trước khi chạy thực nghiệm.")


if __name__ == "__main__":
    main()
