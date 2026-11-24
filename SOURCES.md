# Nguồn dữ liệu và công thức

Bộ chương trình dùng định dạng dữ liệu MATPOWER version 2 và triển khai quy ước DC tương ứng.

1. MATPOWER, `makeBdc`: https://matpower.org/documentation/ref-manual/legacy/functions/makeBdc.html
   - Quy ước `Pbus = Bbus * Va + Pbusinj` và `Pf = Bf * Va + Pfinj`.
   - Hệ số nhánh có xét trạng thái nhánh, điện kháng, tap ratio và phase shift.
2. MATPOWER, Data File Format: https://matpower.app/manual/matpower/DataFileFormat.html
   - Cấu trúc các trường `baseMVA`, `bus`, `gen`, `branch` và ý nghĩa các cột dữ liệu.
3. MATPOWER repository: https://github.com/MATPOWER/matpower/tree/master/data
   - Nguồn các case `case14.m`, `case30.m`, `case57.m`, `case118.m`, `case300.m`.
4. Zimmerman, R. D., Murillo-Sanchez, C. E., & Thomas, R. J. (2011). MATPOWER: Steady-State Operations, Planning and Analysis Tools for Power Systems Research and Education. IEEE Transactions on Power Systems, 26(1), 12–19. DOI: 10.1109/TPWRS.2010.2051168.

Ghi chú mô hình: công suất nút quy định là tổng công suất phát của các máy phát online trừ công suất phụ tải. Nút REF giữ góc pha cố định, phương trình ở nút REF bị loại khỏi hệ rút gọn; công suất nút REF được xác định từ trạng thái cân bằng sau khi giải.
