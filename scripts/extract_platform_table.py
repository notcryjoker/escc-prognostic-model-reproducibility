#!/usr/bin/env python3
"""
流式解压 GPL18109_family.soft.gz，边解边抽平台注释表。
原理: gzip 为顺序流，平台表位于文件前部；即使文件尾部损坏，
只要开头有效即可完整提取 !platform_table_begin ~ !platform_table_end 之间的内容。
"""
import zlib, sys, os

src = sys.argv[1] if len(sys.argv) > 1 else "GPL18109_family.soft.gz"
dst = "GPL18109_platform_table.txt"

d = zlib.decompressobj(16 + zlib.MAX_WBITS)
fh = open(src, "rb")
out = open(dst, "w", encoding="utf-8", errors="replace")

buf = b""
in_table = False
n_lines = 0
n_table = 0
done = False
total_in = 0

while not done:
    chunk = fh.read(1 << 20)
    if not chunk:
        break
    total_in += len(chunk)
    try:
        buf += d.decompress(chunk)
    except zlib.error as e:
        print(f"[WARN] 解压中断于 {total_in/1048576:.1f} MB: {e}", flush=True)
        break
    lines = buf.split(b"\n")
    buf = lines[-1]
    for ln in lines[:-1]:
        s = ln.decode("utf-8", errors="replace")
        n_lines += 1
        if not in_table and s.startswith("!platform_table_begin"):
            in_table = True
            continue
        if in_table:
            if s.startswith("!platform_table_end"):
                out.close()
                print(f"[OK] 平台表提取完成: {n_table} 行, 已消费压缩流 {total_in/1048576:.1f} MB", flush=True)
                done = True
                break
            out.write(s + "\n")
            n_table += 1
            if n_table % 5000 == 0:
                print(f"  ... {n_table} 行", flush=True)

if not done:
    # 冲刷剩余
    for ln in buf.split(b"\n"):
        s = ln.decode("utf-8", errors="replace")
        if in_table and not s.startswith("!platform_table_end"):
            out.write(s + "\n"); n_table += 1
    out.close()
    print(f"[PARTIAL] 流提前结束: 共解析 {n_lines} 行, 平台表 {n_table} 行, "
          f"已消费 {total_in/1048576:.1f} MB", flush=True)

fh.close()
print(f"输出: {dst}  ({os.path.getsize(dst)} bytes)", flush=True)
