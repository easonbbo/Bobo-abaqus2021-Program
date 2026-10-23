#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
从当前打开的 Abaqus ODB 提取环件轧制历史输出节点坐标，计算椭圆度并导出为 CSV。

用法：
  在 Abaqus/CAE 中打开 ODB，然后 File > Run Script 运行本脚本。
  脚本自动获取当前会话中打开的 ODB（只读取，不会关闭它），
  自动识别 8 个圈（SA~SH）每圈 24 个节点的历史输出（COOR1/COOR3），
  逐时间点拟合椭圆，输出每个圈一个 CSV 文件。

输出目录：C:\Users\30557\Desktop\themost commonly uesdfiles\EXCEL
"""

import csv
import os
import re
import sys

import numpy as np

# ================== 可修改配置 ==================
OUTPUT_DIR = r'C:\Users\30557\Desktop\themost commonly uesdfiles\EXCEL'  # 输出目录

# 历史输出变量名（X / Z 坐标）
VAR_X = 'COOR1'
VAR_Z = 'COOR3'

# 8 个圈，每圈 24 个点；节点标签从 ODB 自动识别，无需手改编号
CIRCLE_NAMES = ['SA', 'SB', 'SC', 'SD', 'SE', 'SF', 'SG', 'SH']
NODE_PER_CIRCLE = 24


# ================== 椭圆拟合函数（Fitzgibbon 直接最小二乘） ==================
def fit_ellipse(points):
    """
    输入: points 为 Nx2 的 numpy 数组，每行 (x, z)
    输出: 长半轴 a, 短半轴 b
    """
    # Fitzgibbon 直接法对近圆形点列数值不稳定（特征值符号会翻转），
    # 失败时回退到 PCA 主轴法：对均匀分布在圆周上的 24 个点，
    # PCA 长短轴即为椭圆长短轴的精确估计。
    def pca_fallback(pts):
        center = pts.mean(axis=0)
        cov = np.cov((pts - center).T, bias=1)
        w = np.sort(np.linalg.eigvalsh(cov))[::-1]
        a = np.sqrt(2.0 * max(w[0], 0.0))
        b = np.sqrt(2.0 * max(w[1], 0.0))
        if a < b:
            a, b = b, a
        return a, b

    x = points[:, 0]
    z = points[:, 1]
    D = np.vstack([x*x, x*z, z*z, x, z, np.ones_like(x)]).T
    S = np.dot(D.T, D)
    C = np.zeros((6, 6))
    C[0, 2] = 2
    C[2, 0] = 2
    C[1, 1] = -1
    try:
        eigval, eigvec = np.linalg.eig(np.dot(np.linalg.inv(S), C))
    except np.linalg.LinAlgError:
        S += np.eye(6) * 1e-10
        eigval, eigvec = np.linalg.eig(np.dot(np.linalg.inv(S), C))
    idx = np.where((eigval > 0) & (np.isreal(eigval)))[0]
    if len(idx) == 0:
        return pca_fallback(points)
    a_vec = np.real(eigvec[:, idx[0]])
    A, B, C_, D_, E_, F = a_vec
    denom = B**2 - 4*A*C_
    if abs(denom) < 1e-12:
        return pca_fallback(points)
    xc = (2*C_*D_ - B*E_) / denom
    zc = (2*A*E_ - B*D_) / denom
    F2 = A*xc**2 + B*xc*zc + C_*zc**2 - F
    theta = 0.5 * np.arctan2(B, A - C_)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    a2 = F2 / (A*cos_t**2 + B*sin_t*cos_t + C_*sin_t**2)
    b2 = F2 / (A*sin_t**2 - B*sin_t*cos_t + C_*cos_t**2)
    a = np.sqrt(abs(a2))
    b = np.sqrt(abs(b2))
    if a < b:
        a, b = b, a
    return a, b


# ================== 获取当前 ODB ==================
def get_current_odb():
    """
    从当前 Abaqus/CAE 会话获取已打开的 ODB。
    注意：只读取 session.odbs，绝不调用 odb.close()，不会关闭用户打开的 ODB。
    """
    try:
        from abaqus import session
        odbs = session.odbs
        if odbs:
            # 取最后打开的一个 ODB
            odb = odbs[list(odbs.keys())[-1]]
            print("使用当前打开的 ODB:", odb.name)
            return odb
        else:
            print("错误：当前会话没有打开的 ODB。")
            return None
    except ImportError:
        # 不在 Abaqus/CAE 环境（例如独立 Python）
        return None


# ================== 从区域名识别节点标签 ==================
def parse_region_label(region_name):
    """
    从历史区域名中提取 (实例名, 节点标签)。
    支持 'Node RING-1.6'、'RING-1.6'、'Node ASSEMBLY.6' 等格式。
    """
    m = re.match(r'^(?:Node\s+)?(.+?)\.(\d+)$', region_name.strip())
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def detect_circles(region_names):
    """
    自动识别 8 圈 x 24 点的节点网格。

    历史输出区域按节点命名（如 'Node RING-1.6'），节点标签构成
    8 列 x 24 行的网格：每圈 24 个节点，标签 = 起始标签 + 行号*步长。

    返回 (lookup, meta)：
      lookup[label] = (圈名, 圈内序号 0~23)
      meta 含识别出的起始标签、步长、实例名、各实例节点数
    识别失败时 lookup 为 None。
    """
    labels_by_inst = {}
    for rname in region_names:
        inst, label = parse_region_label(rname)
        if label is None:
            continue
        labels_by_inst.setdefault(inst, set()).add(label)

    meta = {
        'by_inst': {k: len(v) for k, v in labels_by_inst.items()},
        'starts': None,
        'step': None,
        'instance': None,
    }

    expected = len(CIRCLE_NAMES) * NODE_PER_CIRCLE
    candidates = [(inst, labels) for inst, labels in labels_by_inst.items()
                  if len(labels) == expected]
    if not candidates:
        return None, meta

    # 优先选择名字含 RING 的实例
    ring_cands = [c for c in candidates if 'ring' in c[0].lower()]
    inst, labels = ring_cands[0] if ring_cands else candidates[0]

    sorted_labels = sorted(labels)
    starts = sorted_labels[:len(CIRCLE_NAMES)]
    step = sorted_labels[len(CIRCLE_NAMES)] - starts[0]

    lookup = {}
    for c, start in enumerate(starts):
        for k in range(NODE_PER_CIRCLE):
            lab = start + k * step
            if lab in labels:
                lookup[lab] = (CIRCLE_NAMES[c], k)
            else:
                print("警告: 期望节点 %d 不存在于历史输出中" % lab)

    meta.update({'starts': starts, 'step': step, 'instance': inst})
    return lookup, meta


# ================== 主流程 ==================
def main():
    # 获取当前打开的 ODB（不关闭）
    odb = get_current_odb()
    if odb is None:
        print("无法自动获取 ODB。请确认脚本在 Abaqus/CAE 中运行且已打开 ODB。")
        sys.exit(1)

    # 使用最后一个分析步
    step_names = list(odb.steps.keys())
    if not step_names:
        sys.exit("ODB 里没有任何分析步")
    last_step_name = step_names[-1]
    step = odb.steps[last_step_name]
    print("使用分析步:", last_step_name)

    # 自动识别 8 圈 x 24 点
    region_names = list(step.historyRegions.keys())
    lookup, meta = detect_circles(region_names)
    if lookup is None:
        print("错误: 未能识别出 %d 个圈 x 每圈 %d 个点的历史输出网格。"
              % (len(CIRCLE_NAMES), NODE_PER_CIRCLE))
        print("各实例的历史节点区域数量:", meta['by_inst'])
        print("前 30 个历史区域名:", region_names[:30])
        sys.exit(1)

    print("识别实例:", meta['instance'])
    print("节点间隔:", meta['step'])
    print("圈 -> 起始节点标签:", dict(zip(CIRCLE_NAMES, meta['starts'])))

    # 数据存储结构
    data = {prefix: {} for prefix in CIRCLE_NAMES}

    # 读取历史输出数据
    for region_name, region in step.historyRegions.items():
        inst, label = parse_region_label(region_name)
        # 只取识别出的那个实例，避免不同实例节点标签相同导致串数据
        if inst != meta['instance'] or label not in lookup:
            continue
        prefix, index = lookup[label]

        if VAR_X not in region.historyOutputs or VAR_Z not in region.historyOutputs:
            print("警告: 节点 %d 缺少 %s 或 %s" % (label, VAR_X, VAR_Z))
            continue

        coord1_data = region.historyOutputs[VAR_X].data
        coord3_data = region.historyOutputs[VAR_Z].data

        times = [t for t, v in coord1_data]
        x_vals = [v for t, v in coord1_data]
        z_vals = [v for t, v in coord3_data]

        data[prefix][index] = {
            'times': times,
            'x': x_vals,
            'z': z_vals,
        }

    # 检查每个圈缺失的节点
    for prefix in CIRCLE_NAMES:
        missing = [i + 1 for i in range(NODE_PER_CIRCLE) if i not in data[prefix]]
        if missing:
            print("警告: 圈 %s 缺少节点: %s" % (prefix, missing))

    # 获取全局时间序列：取时间点最多的节点，避免个别节点帧数少导致数据截断
    global_times = None
    for prefix in CIRCLE_NAMES:
        for idx in data[prefix]:
            cand_times = data[prefix][idx]['times']
            if global_times is None or len(cand_times) > len(global_times):
                global_times = cand_times
    if global_times is None:
        sys.exit("未能获取任何时间数据")

    print("共 %d 个时间点" % len(global_times))

    # 创建输出目录（Python 2 的 makedirs 没有 exist_ok 参数，先判断再创建）
    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # 为每个圈生成 CSV 文件
    for prefix in CIRCLE_NAMES:
        csv_path = os.path.join(OUTPUT_DIR, '%s.csv' % prefix)
        # Python 2 写二进制文件，Python 3 写 UTF-8 文本文件
        if sys.version_info[0] >= 3:
            csv_file = open(csv_path, 'w', newline='', encoding='utf-8')
        else:
            csv_file = open(csv_path, 'wb')
        with csv_file as f:
            writer = csv.writer(f)
            # 表头
            header = ['t']
            for i in range(1, NODE_PER_CIRCLE + 1):
                header.append('%s%d_X' % (prefix, i))
                header.append('%s%d_Z' % (prefix, i))
            header += ['MAX-Diameter', 'MIN-Diameter', 'ovality']
            writer.writerow(header)

            # 逐时间点处理
            for t_idx, t in enumerate(global_times):
                points = []
                row = [t]
                valid = True
                for i in range(NODE_PER_CIRCLE):
                    # 防止个别节点时间点数少于全局时间序列导致索引越界
                    if (i in data[prefix]
                            and t_idx < len(data[prefix][i]['x'])
                            and t_idx < len(data[prefix][i]['z'])):
                        x = data[prefix][i]['x'][t_idx]
                        z = data[prefix][i]['z'][t_idx]
                        points.append([x, z])
                        row.append(x)
                        row.append(z)
                    else:
                        points.append([np.nan, np.nan])
                        row.append(None)
                        row.append(None)
                        valid = False
                if not valid or len(points) < 6:
                    row += [None, None, None]
                    writer.writerow(row)
                    continue

                points = np.array(points)
                mask = ~np.isnan(points).any(axis=1)
                points_clean = points[mask]
                if len(points_clean) < 6:
                    row += [None, None, None]
                    writer.writerow(row)
                    continue

                try:
                    a, b = fit_ellipse(points_clean)
                    Dmax = 2 * a
                    Dmin = 2 * b
                    oval = Dmax - Dmin
                    row += [Dmax, Dmin, oval]
                except Exception as e:
                    print("警告: 时间 %f 圈 %s 拟合失败: %s" % (t, prefix, str(e)))
                    row += [None, None, None]
                writer.writerow(row)

        print("已生成:", csv_path)

    print("完成！所有 CSV 文件已保存至:", OUTPUT_DIR)


if __name__ == '__main__':
    main()
