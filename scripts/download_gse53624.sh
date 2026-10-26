#!/bin/bash
# GSE53624 系列矩阵多段并行下载（零填充命名，避免字典序拼接错误）
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data/geo"

URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53624/matrix/GSE53624_series_matrix.txt.gz"
SIZE=71937810
HAVE=26800000          # 已下载部分（小于实际，脚本会重新校验）
N=16

# 若已有部分文件，先记录长度以便复用
if [ -f GSE53624_series_matrix.txt.gz ]; then
  HAVE=$(stat -c%s GSE53624_series_matrix.txt.gz)
fi
echo "已有 $HAVE / $SIZE 字节"

# 清理旧的零填充分片
rm -f gse_part_*.bin

# 剩余部分分段下载
REM=$((SIZE - HAVE))
CHUNK=$(( (REM + N - 1) / N ))
echo "剩余 $REM 字节, 分 $N 段, 每段 $CHUNK 字节"

for i in $(seq 0 $((N-1))); do
  printf -v IDX "%02d" $i
  START=$((HAVE + i*CHUNK))
  END=$((HAVE + (i+1)*CHUNK - 1))
  [ $END -ge $SIZE ] && END=$((SIZE-1))
  [ $START -gt $END ] && continue
  curl -s -r ${START}-${END} --max-time 2400 --retry 8 --retry-all-errors \
       -o gse_part_${IDX}.bin "$URL" &
done
wait

# 按数字序拼接（零填充保证字典序 == 数字序）
cat gse_part_*.bin >> GSE53624_series_matrix.txt.gz
rm -f gse_part_*.bin

CUR=$(stat -c%s GSE53624_series_matrix.txt.gz)
echo "拼接后大小: $CUR / $SIZE"
if gzip -t GSE53624_series_matrix.txt.gz 2>/dev/null; then
  echo "GSE53624_OK $CUR"
else
  echo "GSE53624_FAIL $CUR"
fi
