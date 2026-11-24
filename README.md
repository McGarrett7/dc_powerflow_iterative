# Thực nghiệm DC Power Flow bằng Jacobi, Gauss–Seidel, SOR và SSOR

Bộ chương trình Python này xây dựng hệ DC Power Flow từ các file case định dạng MATPOWER, lưu ma trận dưới dạng CSR, chạy bốn phương pháp lặp và xuất CSV/biểu đồ để phục vụ Chương 5 của tiểu luận.

## 1. Cấu trúc thư mục

```text
dc_powerflow_iterative/
├── data/
│   └── case14.m                 # IEEE 14-bus, đã kèm sẵn
├── src/
│   ├── matpower_parser.py       # Đọc file MATPOWER dạng .m
│   ├── dc_pf.py                 # Xây dựng Bbus, B', P và công suất nhánh
│   └── iterative_solvers.py     # Jacobi, Gauss-Seidel, SOR, SSOR
├── tests/
│   └── test_solvers.py
├── results/                     # Kết quả CSV và PNG
├── download_cases.py            # Tải case30/57/118/300 từ MATPOWER
├── run_experiments.py            # Chạy thực nghiệm
└── requirements.txt
```

## 2. Cài đặt môi trường

Khuyến nghị Python 3.10 trở lên.

```bash
python -m venv .venv
```

Windows:

```powershell
.venv\Scripts\Activate.ps1
```

Ubuntu/Linux:

```bash
source .venv/bin/activate
```

Cài thư viện:

```bash
python -m pip install -r requirements.txt
```

## 3. Dữ liệu IEEE bus

File `case14.m` được kèm sẵn trong thư mục `data/`. Để tải thêm các case IEEE 30, 57, 118 và 300 bus từ repository MATPOWER chính thức, chạy:

```bash
python download_cases.py
```

Các file tải về là `data/case30.m`, `data/case57.m`, `data/case118.m`, `data/case300.m`.

Nếu máy không truy cập Internet, tải thủ công các file cùng tên từ:

`https://github.com/MATPOWER/matpower/tree/master/data`

và đặt vào thư mục `data/`.

## 4. Chạy thử IEEE 14-bus

```bash
python run_experiments.py --cases case14 --repeats 3 --tol 1e-8
```

## 5. Chạy tất cả các hệ thống IEEE được cấu hình

Sau khi tải đủ dữ liệu:

```bash
python run_experiments.py --cases case14 case30 case57 case118 case300 --repeats 3 --tol 1e-8
```

Tham số mặc định:

- Ngưỡng phần dư tương đối: `1e-8`.
- Số vòng lặp tối đa: `200000`.
- Số lần đo thời gian: `3`; báo cáo trung vị thời gian.
- Khảo sát SOR/SSOR: `omega` từ `0.80` đến `1.95`, bước `0.05`.
- Tiêu chuẩn chọn `omega`: số vòng lặp nhỏ nhất trong các giá trị khảo sát hội tụ; nếu bằng nhau, chọn thời gian đo thấp hơn.

Có thể thay đổi tham số, ví dụ:

```bash
python run_experiments.py --cases case14 case118 --tol 1e-10 --max-iter 300000 --repeats 5 --omega-step 0.02
```

Để bỏ khảo sát omega (chương trình dùng omega=1.0 cho SOR/SSOR):

```bash
python run_experiments.py --cases case14 --no-omega-sweep
```

## 6. Các file kết quả

Trong `results/`:

- `matrix_statistics.csv`: số nút, số góc pha ẩn, số nhánh, `nnz`, mật độ và dung lượng CSR so với ma trận đặc.
- `solver_summary.csv`: số vòng lặp, phần dư cuối, thời gian trung vị và sai số góc pha so với nghiệm tham chiếu.
- `omega_sweep.csv`: kết quả khảo sát từng giá trị omega của SOR và SSOR.
- `bus_angles_case*.csv`: góc pha và công suất nút.
- `branch_flows_case*.csv`: công suất tác dụng ở đầu from-to và to-from của từng nhánh.
- `convergence_history_case*.csv`: lịch sử phần dư của các phương pháp.
- `iterations_by_case.png`, `runtime_by_case.png`: biểu đồ so sánh tổng hợp.
- `omega_sweep_case*.png`: biểu đồ ảnh hưởng omega.
- `convergence_case*.png`: đường cong hội tụ.
- `experiment_config.json`: tham số và môi trường thực nghiệm.

## 7. Mô hình toán học và quy ước

Chương trình áp dụng quy ước DC của MATPOWER:

`Pbus = Bbus @ Va + Pbusinj`

`Pbranch_from = Bf @ Va + Pfinj`

Trong đó góc `Va` tính bằng radian, công suất nội bộ tính theo p.u. Khi xuất kết quả nhánh, công suất được nhân với `baseMVA` để chuyển sang MW.

- Công suất nút quy định: tổng `Pg` của các máy phát đang hoạt động tại nút trừ `Pd`.
- Nút REF được giữ cố định góc pha theo dữ liệu case; hàng/cột nút REF được loại khỏi hệ cần giải.
- Có xét trạng thái nhánh, tỷ số tap khác 0 và góc lệch pha máy biến áp.
- Không dùng giới hạn công suất nhánh hoặc giới hạn điện áp trong DC Power Flow tuyến tính cơ bản.
- Nghiệm tham chiếu được tính bằng sparse direct solver (`scipy.sparse.linalg.spsolve`) chỉ để kiểm chứng độ chính xác; bốn kết quả chính là từ các phương pháp lặp.

## 8. Kiểm thử

```bash
python -m pytest tests
```

Test kiểm tra bốn thuật toán trên một ma trận SPD nhỏ và kiểm tra bộ dựng mô hình DC cho IEEE 14-bus.

## 9. Lưu ý khi đưa kết quả vào báo cáo

1. Không điền số liệu từ ví dụ hoặc từ case khác vào bảng IEEE bus.
2. Dùng cùng tolerance, `x0`, cấu hình máy và số lần lặp đo thời gian cho tất cả phương pháp.
3. Báo cáo rõ SOR/SSOR dùng omega nào và omega được chọn theo tiêu chí nào.
4. Trình bày riêng thời gian xây dựng ma trận và thời gian giải hệ nếu đo cả hai.
5. DC Power Flow là mô hình tuyến tính xấp xỉ; công suất nhánh không bao gồm tổn thất công suất tác dụng.
6. Các case MATPOWER là bộ dữ liệu kiểm thử chuẩn, không đại diện cho trạng thái vận hành hiện tại của một lưới điện thực tế.

## 10. Nguồn dữ liệu và tài liệu kỹ thuật

- MATPOWER repository: https://github.com/MATPOWER/matpower
- MATPOWER `makeBdc`: https://matpower.org/documentation/ref-manual/legacy/functions/makeBdc.html
- MATPOWER case format: https://matpower.org/documentation/ref-manual/legacy/functions/caseformat.html
- MATPOWER paper: Zimmerman, Murillo-Sanchez & Thomas (2011), IEEE Transactions on Power Systems, 26(1), 12–19, DOI: 10.1109/TPWRS.2010.2051168.
