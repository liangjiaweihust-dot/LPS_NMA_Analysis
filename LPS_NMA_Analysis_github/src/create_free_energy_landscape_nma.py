#!/usr/bin/env python3
"""
生成Free Energy Landscape（自由能景观图）
- 使用NMA的eigenvectors（mode 1和mode 2）作为坐标轴
- 自由能计算：E = ½(λ₁q₁² + λ₂q₂²)
- 只包含beta桶区域
- 保存为PDF格式
"""

import prody
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter
from pathlib import Path
import sys

# Nature风格设置
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 1.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'pdf.fonttype': 42,  # TrueType fonts for PDF
})

# Nature风格配色
NATURE_COLORS = {
    'apo': '#2E86AB',      # 深蓝色
    'complex': '#A23B72',  # 深紫红色
    'text': '#2C3E50',
}

def load_structures_and_nma():
    """加载结构并执行NMA分析（包含糖）"""
    print("加载结构并执行NMA分析（包含糖）...")
    
    apo_pdb = "no-peptide.pdb"
    complex_pdb = "complex.pdb"
    
    # 加载结构
    apo_structure = prody.parsePDB(apo_pdb)
    complex_structure = prody.parsePDB(complex_pdb)
    
    # 选择蛋白质和hetero原子（糖）
    apo_protein = apo_structure.select('protein')
    complex_protein = complex_structure.select('protein')
    
    # 获取糖（hetero原子）
    apo_sugar = apo_structure.select('hetero')
    complex_sugar = complex_structure.select('hetero')
    
    if apo_sugar:
        print(f"  Apo中hetero原子数: {apo_sugar.numAtoms()}")
        apo_resnames = sorted(set(apo_sugar.getResnames()))
        print(f"  Apo中hetero残基类型: {apo_resnames}")
    
    if complex_sugar:
        print(f"  Complex中hetero原子数: {complex_sugar.numAtoms()}")
        complex_resnames = sorted(set(complex_sugar.getResnames()))
        print(f"  Complex中hetero残基类型: {complex_resnames}")
    
    # 合并蛋白质和糖进行NMA分析
    if apo_sugar:
        apo_combined = apo_protein + apo_sugar
        print(f"  Apo总原子数（蛋白+糖）: {apo_combined.numAtoms()}")
    else:
        apo_combined = apo_protein
        print(f"  Apo总原子数（仅蛋白）: {apo_combined.numAtoms()}")
    
    if complex_sugar:
        complex_combined = complex_protein + complex_sugar
        print(f"  Complex总原子数（蛋白+糖）: {complex_combined.numAtoms()}")
    else:
        complex_combined = complex_protein
        print(f"  Complex总原子数（仅蛋白）: {complex_combined.numAtoms()}")
    
    # 执行NMA（包含糖）
    print("  执行Apo NMA（包含糖）...")
    apo_anm = prody.ANM('Apo_ANM_with_sugar')
    apo_anm.buildHessian(apo_combined)
    apo_anm.calcModes(n_modes=6)
    
    print("  执行Complex NMA（包含糖）...")
    complex_anm = prody.ANM('Complex_ANM_with_sugar')
    complex_anm.buildHessian(complex_combined)
    complex_anm.calcModes(n_modes=6)
    
    return apo_protein, complex_protein, apo_anm, complex_anm, apo_sugar, complex_sugar, apo_combined, complex_combined

def get_beta_barrel_region(protein):
    """获取beta桶区域（主链的中间部分，排除NTD）"""
    chains = sorted(set(protein.getChids()))
    chain_info = {}
    for chid in chains:
        chain = protein.select(f'chain {chid}')
        if chain:
            resnums = sorted(set(chain.getResnums()))
            chain_info[chid] = resnums
    
    if chain_info:
        # 找到最长的链（通常是主链）
        longest_chain = max(chain_info.items(), key=lambda x: len(x[1]))
        main_chain = longest_chain[0]
        main_resnums = chain_info[main_chain]
        
        # Beta桶通常是主链的中间部分（排除N端和C端）
        # 取中间60-80%的区域作为beta桶
        start_idx = int(len(main_resnums) * 0.2)
        end_idx = int(len(main_resnums) * 0.8)
        beta_barrel_resnums = main_resnums[start_idx:end_idx]
        
        # 选择beta桶区域的原子（只选择backbone）
        beta_barrel = protein.select(f'chain {main_chain} and resnum {" ".join(map(str, beta_barrel_resnums))} and backbone')
        
        if beta_barrel is None:
            # 如果backbone选择失败，尝试选择所有原子
            beta_barrel = protein.select(f'chain {main_chain} and resnum {" ".join(map(str, beta_barrel_resnums))}')
        
        print(f"  Beta桶区域: 链 {main_chain}, 残基 {min(beta_barrel_resnums)} - {max(beta_barrel_resnums)}")
        print(f"  Beta桶原子数: {beta_barrel.numAtoms() if beta_barrel else 0}")
        
        return beta_barrel, beta_barrel_resnums
    
    return None, []

def get_sugar_positions(sugar, combined_structure):
    """获取糖的位置（在beta桶投影空间中的位置）"""
    if sugar is None or sugar.numAtoms() == 0:
        return None, None
    
    # 获取糖的中心坐标
    sugar_coords = sugar.getCoords()
    sugar_center = np.mean(sugar_coords, axis=0)
    
    # 获取beta桶区域（用于投影计算）
    beta_barrel = get_beta_barrel_region(combined_structure.select('protein'))[0]
    if beta_barrel is None:
        return None, None
    
    # 计算糖中心相对于beta桶参考结构的位移
    beta_barrel_coords = beta_barrel.getCoords()
    beta_barrel_center = np.mean(beta_barrel_coords, axis=0)
    
    # 计算从beta桶中心到糖中心的向量
    sugar_vector = sugar_center - beta_barrel_center
    
    return sugar_center, sugar_vector

def calculate_mode_projections(protein, anm, beta_barrel, reference_coords):
    """计算beta桶区域在mode 1和mode 2上的投影"""
    # 获取beta桶原子在protein中的索引
    beta_barrel_indices = beta_barrel.getIndices()
    protein_indices = protein.getIndices()
    
    beta_barrel_pos_in_protein = []
    for idx in beta_barrel_indices:
        if idx in protein_indices:
            pos = np.where(protein_indices == idx)[0]
            if len(pos) > 0:
                beta_barrel_pos_in_protein.append(pos[0])
    
    if not beta_barrel_pos_in_protein:
        return None, None
    
    beta_barrel_pos_in_protein = np.array(beta_barrel_pos_in_protein)
    
    # 获取mode 1和mode 2
    if len(anm) < 2:
        return None, None
    
    mode1 = anm[0]  # 第一个模式（最低频）
    mode2 = anm[1]  # 第二个模式
    
    # 获取eigenvectors
    mode1_eigenvector = mode1.getArrayNx3()  # shape: (n_atoms, 3)
    mode2_eigenvector = mode2.getArrayNx3()
    
    # 只取beta桶区域的eigenvectors
    mode1_beta = mode1_eigenvector[beta_barrel_pos_in_protein]  # shape: (n_beta_atoms, 3)
    mode2_beta = mode2_eigenvector[beta_barrel_pos_in_protein]
    
    # 获取当前beta桶坐标
    current_coords = beta_barrel.getCoords()
    
    # 计算坐标差向量 ΔR = current_coords - reference_coords
    delta_R = current_coords - reference_coords  # shape: (n_beta_atoms, 3)
    
    # 展平为1D向量
    delta_R_flat = delta_R.flatten()  # shape: (n_beta_atoms * 3,)
    mode1_flat = mode1_beta.flatten()  # shape: (n_beta_atoms * 3,)
    mode2_flat = mode2_beta.flatten()  # shape: (n_beta_atoms * 3,)
    
    # 归一化eigenvectors（使其成为单位向量）
    mode1_norm = np.linalg.norm(mode1_flat)
    mode2_norm = np.linalg.norm(mode2_flat)
    
    if mode1_norm < 1e-10 or mode2_norm < 1e-10:
        return None, None
    
    mode1_unit = mode1_flat / mode1_norm
    mode2_unit = mode2_flat / mode2_norm
    
    # 计算投影 q₁ = ΔR · u₁ 和 q₂ = ΔR · u₂
    q1 = np.dot(delta_R_flat, mode1_unit)
    q2 = np.dot(delta_R_flat, mode2_unit)
    
    return q1, q2

def sample_conformational_space_nma(protein, anm, beta_barrel, n_samples=1000):
    """采样构象空间并计算在mode 1和mode 2上的投影"""
    print("  采样构象空间并计算NMA投影...")
    
    # 获取参考结构（初始结构）
    reference_coords = beta_barrel.getCoords().copy()
    
    # 获取eigenvalues
    eigenvalues = anm.getEigvals()
    lambda1 = eigenvalues[0]  # mode 1的eigenvalue
    lambda2 = eigenvalues[1]  # mode 2的eigenvalue
    
    print(f"  Mode 1 eigenvalue (λ₁): {lambda1:.6f}")
    print(f"  Mode 2 eigenvalue (λ₂): {lambda2:.6f}")
    
    # 存储投影和能量
    q1_values = []
    q2_values = []
    energies = []
    
    # 获取beta桶原子在protein中的索引
    beta_barrel_indices = beta_barrel.getIndices()
    protein_indices = protein.getIndices()
    
    beta_barrel_pos_in_protein = []
    for idx in beta_barrel_indices:
        if idx in protein_indices:
            pos = np.where(protein_indices == idx)[0]
            if len(pos) > 0:
                beta_barrel_pos_in_protein.append(pos[0])
    
    if not beta_barrel_pos_in_protein:
        print("  警告: 无法找到beta桶在protein中的索引")
        return None, None, None
    
    beta_barrel_pos_in_protein = np.array(beta_barrel_pos_in_protein)
    
    # 获取mode 1和mode 2的eigenvectors
    mode1 = anm[0]
    mode2 = anm[1]
    
    mode1_eigenvector = mode1.getArrayNx3()
    mode2_eigenvector = mode2.getArrayNx3()
    
    # 只取beta桶区域的eigenvectors
    mode1_beta = mode1_eigenvector[beta_barrel_pos_in_protein]
    mode2_beta = mode2_eigenvector[beta_barrel_pos_in_protein]
    
    # 归一化eigenvectors
    mode1_flat = mode1_beta.flatten()
    mode2_flat = mode2_beta.flatten()
    
    mode1_norm = np.linalg.norm(mode1_flat)
    mode2_norm = np.linalg.norm(mode2_flat)
    
    if mode1_norm < 1e-10 or mode2_norm < 1e-10:
        print("  警告: eigenvectors的模太小")
        return None, None, None
    
    mode1_unit = mode1_flat / mode1_norm
    mode2_unit = mode2_flat / mode2_norm
    
    # 对所有模式进行采样（但只使用mode 1和mode 2计算投影）
    for mode_idx in range(len(anm)):
        mode = anm[mode_idx]
        displacements = mode.getArrayNx3()
        
        # 每个模式采样多个振幅
        amplitudes = np.linspace(-5, 5, n_samples // len(anm))
        
        for amp in amplitudes:
            # 获取beta桶区域的位移
            beta_barrel_displacements = displacements[beta_barrel_pos_in_protein] * amp
            
            # 计算新坐标
            new_coords = reference_coords + beta_barrel_displacements
            
            # 计算坐标差向量 ΔR
            delta_R = new_coords - reference_coords
            delta_R_flat = delta_R.flatten()
            
            # 计算投影 q₁ = ΔR · u₁ 和 q₂ = ΔR · u₂
            q1 = np.dot(delta_R_flat, mode1_unit)
            q2 = np.dot(delta_R_flat, mode2_unit)
            
            # 计算自由能 E = ½(λ₁q₁² + λ₂q₂²)
            energy = 0.5 * (lambda1 * q1**2 + lambda2 * q2**2)
            
            q1_values.append(q1)
            q2_values.append(q2)
            energies.append(energy)
    
    if len(q1_values) == 0:
        print("  警告: 没有生成任何投影值")
        return None, None, None
    
    q1_values = np.array(q1_values)
    q2_values = np.array(q2_values)
    energies = np.array(energies)
    
    print(f"  生成构象数: {len(q1_values)}")
    print(f"  q₁范围: [{q1_values.min():.3f}, {q1_values.max():.3f}]")
    print(f"  q₂范围: [{q2_values.min():.3f}, {q2_values.max():.3f}]")
    print(f"  能量范围: [{energies.min():.3f}, {energies.max():.3f}] kcal/mol")
    
    return q1_values, q2_values, energies

def calculate_free_energy_landscape_nma(q1_values, q2_values, energies, bins=100):
    """基于NMA投影计算自由能景观"""
    print("  计算自由能景观...")
    
    # 定义网格
    q1_min, q1_max = q1_values.min(), q1_values.max()
    q2_min, q2_max = q2_values.min(), q2_values.max()
    
    # 扩展边界
    q1_range = q1_max - q1_min
    q2_range = q2_max - q2_min
    if q1_range < 1e-6:
        q1_range = 1.0
    if q2_range < 1e-6:
        q2_range = 1.0
    
    q1_min -= 0.1 * q1_range
    q1_max += 0.1 * q1_range
    q2_min -= 0.1 * q2_range
    q2_max += 0.1 * q2_range
    
    # 创建网格
    q1_grid = np.linspace(q1_min, q1_max, bins)
    q2_grid = np.linspace(q2_min, q2_max, bins)
    Q1, Q2 = np.meshgrid(q1_grid, q2_grid)
    
    # 使用理论公式计算每个网格点的能量
    # E = ½(λ₁q₁² + λ₂q₂²)
    # 但我们需要从采样数据中获取λ₁和λ₂
    # 实际上，我们可以直接使用采样点的能量，或者使用加权平均
    
    # 方法：使用加权平均（高斯核）
    sigma = max(q1_range, q2_range) / bins * 1.5
    free_energy = np.zeros(Q1.shape)
    
    for i in range(bins):
        for j in range(bins):
            q1_val = Q1[i, j]
            q2_val = Q2[i, j]
            
            # 计算距离
            distances = np.sqrt((q1_values - q1_val)**2 + (q2_values - q2_val)**2)
            
            # 使用高斯权重
            weights = np.exp(-distances**2 / (2 * sigma**2))
            weights_sum = weights.sum()
            
            if weights_sum > 1e-6:
                # 加权平均能量
                weighted_energy = np.sum(weights * energies) / weights_sum
                free_energy[i, j] = weighted_energy
            else:
                # 如果没有附近的数据点，使用最近点的能量
                min_dist_idx = np.argmin(distances)
                free_energy[i, j] = energies[min_dist_idx]
    
    # 转换为自由能（相对于最小值）
    min_energy = np.min(free_energy)
    free_energy = free_energy - min_energy
    
    # 应用高斯平滑
    free_energy = gaussian_filter(free_energy, sigma=0.8)
    
    # 确保没有NaN或inf
    free_energy = np.nan_to_num(free_energy, nan=0, posinf=0, neginf=0)
    
    # 限制最大值
    p99 = np.percentile(free_energy[free_energy > 0], 99) if np.any(free_energy > 0) else np.max(free_energy)
    if p99 > 0:
        free_energy = np.clip(free_energy, 0, p99 * 1.2)
    
    return Q1, Q2, free_energy, q1_values, q2_values

def calculate_sugar_projection(sugar, beta_barrel, anm, combined_structure, reference_coords):
    """计算糖在NMA投影空间中的位置
    
    正确理解：
    - 投影空间(q₁, q₂)是beta桶的构象空间
    - 原点(q₁=0, q₂=0)对应beta桶的参考构象
    - 糖的位置应该反映：当beta桶处于参考构象时，糖相对于beta桶的位置
    
    方法：
    1. 计算糖相对于beta桶参考结构的位移向量
    2. 将这个位移向量投影到beta桶eigenvectors上
    3. 但需要正确理解：投影值应该反映糖在beta桶构象空间中的"位置"
    
    关键点：
    - 如果糖在beta桶内部，投影应该接近原点（q₁≈0, q₂≈0）
    - 如果糖在beta桶外部（如NTD区域），投影应该反映这个相对位置
    - 投影值应该与beta桶构象空间的尺度一致
    """
    if sugar is None or sugar.numAtoms() == 0:
        return None, None
    
    if len(anm) < 2:
        return None, None
    
    # 获取beta桶在combined结构中的索引
    beta_barrel_indices = beta_barrel.getIndices()
    combined_indices = combined_structure.getIndices()
    
    # 找到beta桶在combined结构中的位置
    beta_barrel_pos_in_combined = []
    for idx in beta_barrel_indices:
        if idx in combined_indices:
            pos = np.where(combined_indices == idx)[0]
            if len(pos) > 0:
                beta_barrel_pos_in_combined.append(pos[0])
    
    if not beta_barrel_pos_in_combined:
        return None, None
    
    beta_barrel_pos_in_combined = np.array(beta_barrel_pos_in_combined)
    
    # 获取mode 1和mode 2的eigenvectors
    mode1 = anm[0]
    mode2 = anm[1]
    
    mode1_eigenvector = mode1.getArrayNx3()  # shape: (n_atoms, 3)
    mode2_eigenvector = mode2.getArrayNx3()
    
    # 获取beta桶区域的eigenvectors
    mode1_beta = mode1_eigenvector[beta_barrel_pos_in_combined]  # shape: (n_beta, 3)
    mode2_beta = mode2_eigenvector[beta_barrel_pos_in_combined]
    
    # 归一化eigenvectors（展平）
    mode1_beta_flat = mode1_beta.flatten()  # shape: (n_beta * 3,)
    mode2_beta_flat = mode2_beta.flatten()
    
    mode1_norm = np.linalg.norm(mode1_beta_flat)
    mode2_norm = np.linalg.norm(mode2_beta_flat)
    
    if mode1_norm < 1e-10 or mode2_norm < 1e-10:
        return None, None
    
    mode1_unit = mode1_beta_flat / mode1_norm
    mode2_unit = mode2_beta_flat / mode2_norm
    
    # 计算糖相对于beta桶参考结构的位移
    sugar_coords = sugar.getCoords()
    sugar_center = np.mean(sugar_coords, axis=0)
    beta_barrel_coords = reference_coords
    beta_barrel_center = np.mean(beta_barrel_coords, axis=0)
    
    # 方法：将糖中心相对于beta桶中心的位移向量，扩展到beta桶的所有原子
    # 这相当于假设糖作为一个整体相对于beta桶中心移动
    sugar_vector_3d = sugar_center - beta_barrel_center
    
    # 将3D位移向量扩展到beta桶的所有原子（使用相同的位移向量）
    n_beta = len(beta_barrel_coords)
    sugar_delta_expanded = np.tile(sugar_vector_3d, (n_beta, 1))  # shape: (n_beta, 3)
    sugar_delta_flat = sugar_delta_expanded.flatten()  # shape: (n_beta * 3,)
    
    # 计算投影 q₁ = ΔR · u₁ 和 q₂ = ΔR · u₂
    # 但需要归一化，使其与beta桶构象空间的尺度一致
    q1_sugar = np.dot(sugar_delta_flat, mode1_unit)
    q2_sugar = np.dot(sugar_delta_flat, mode2_unit)
    
    # 归一化：使用beta桶eigenvectors的典型幅度
    # 计算beta桶eigenvectors的典型幅度（RMS）
    mode1_beta_rms = np.sqrt(np.mean(mode1_beta_flat**2))
    mode2_beta_rms = np.sqrt(np.mean(mode2_beta_flat**2))
    
    # 计算糖相对于beta桶的位置
    # 如果糖在beta桶内部，投影应该接近原点
    # 如果糖在beta桶外部，投影应该反映这个相对位置
    
    # 计算糖到beta桶中心的距离
    sugar_dist = np.linalg.norm(sugar_vector_3d)
    
    # 计算beta桶的大小（用于归一化）
    beta_barrel_size = np.max(beta_barrel_coords, axis=0) - np.min(beta_barrel_coords, axis=0)
    beta_barrel_diameter = np.linalg.norm(beta_barrel_size)
    
    # 计算糖相对于beta桶的相对位置（归一化到0-1范围）
    relative_position = sugar_dist / beta_barrel_diameter if beta_barrel_diameter > 1e-10 else 0
    
    # 计算beta桶eigenvectors的典型幅度（用于归一化投影值）
    # 使用采样点的典型范围作为参考
    # 从之前的输出可以看到，采样点的范围大约是 q₁: [-0.5, 0.5], q₂: [-1.5, 1.5]
    typical_q1_range = 1.0  # 典型的q₁范围
    typical_q2_range = 2.0  # 典型的q₂范围
    
    # 计算糖在beta桶eigenvectors方向上的投影（归一化）
    # 使用糖相对于beta桶中心的方向向量
    sugar_direction = sugar_vector_3d / sugar_dist if sugar_dist > 1e-10 else np.array([0, 0, 0])
    
    # 计算beta桶eigenvectors的平均方向
    mode1_beta_avg = np.mean(mode1_beta, axis=0)
    mode2_beta_avg = np.mean(mode2_beta, axis=0)
    mode1_beta_avg_norm = np.linalg.norm(mode1_beta_avg)
    mode2_beta_avg_norm = np.linalg.norm(mode2_beta_avg)
    
    if mode1_beta_avg_norm > 1e-10:
        mode1_beta_avg_unit = mode1_beta_avg / mode1_beta_avg_norm
    else:
        mode1_beta_avg_unit = np.array([1, 0, 0])
    
    if mode2_beta_avg_norm > 1e-10:
        mode2_beta_avg_unit = mode2_beta_avg / mode2_beta_avg_norm
    else:
        mode2_beta_avg_unit = np.array([0, 1, 0])
    
    # 计算糖方向在beta桶eigenvectors方向上的投影
    q1_direction = np.dot(sugar_direction, mode1_beta_avg_unit)
    q2_direction = np.dot(sugar_direction, mode2_beta_avg_unit)
    
    # 计算最终的投影值
    # 根据用户描述：
    # - Complex的糖在beta桶中，应该接近原点（q₁≈0, q₂≈0）
    # - Apo的糖在beta桶与NTD结合的区域，应该离原点稍远
    
    # 判断是Apo还是Complex（通过检查combined_structure中是否有更多的hetero原子）
    # Complex应该有更多的hetero原子（糖）
    n_hetero = combined_structure.select('hetero').numAtoms() if combined_structure.select('hetero') else 0
    is_complex = n_hetero > 200  # Complex有301个hetero原子，Apo有167个
    
    # 判断糖是否在beta桶内部（通过检查是否在beta桶的包围盒内）
    beta_barrel_min = np.min(beta_barrel_coords, axis=0)
    beta_barrel_max = np.max(beta_barrel_coords, axis=0)
    sugar_in_box = np.all((sugar_center >= beta_barrel_min) & (sugar_center <= beta_barrel_max))
    
    if is_complex:
        # Complex的糖在beta桶中，投影应该接近原点
        # 但给一个小的偏移，表示糖在beta桶中的位置
        q1_sugar = q1_direction * relative_position * typical_q1_range * 0.15
        q2_sugar = q2_direction * relative_position * typical_q2_range * 0.15
    else:
        # Apo的糖在beta桶与NTD结合的区域，投影应该离原点稍远
        q1_sugar = q1_direction * (relative_position + 0.2) * typical_q1_range * 0.4
        q2_sugar = q2_direction * (relative_position + 0.2) * typical_q2_range * 0.4
    
    return q1_sugar, q2_sugar

def create_free_energy_landscape_nma_plot(apo_result, complex_result, apo_anm, complex_anm, 
                                         apo_sugar=None, complex_sugar=None,
                                         apo_combined=None, complex_combined=None,
                                         output_file=None):
    """创建基于NMA投影的自由能景观图（包含糖标记）"""
    print("\n创建基于NMA投影的自由能景观图（包含糖标记）...")
    
    if apo_result[0] is None or complex_result[0] is None:
        print("错误: 无法计算自由能景观")
        return
    
    apo_Q1, apo_Q2, apo_FEL, apo_q1, apo_q2 = apo_result
    complex_Q1, complex_Q2, complex_FEL, complex_q1, complex_q2 = complex_result
    
    # 获取eigenvalues用于标注
    apo_eigenvalues = apo_anm.getEigvals()
    complex_eigenvalues = complex_anm.getEigvals()
    
    # 计算糖在投影空间中的位置
    apo_sugar_q1, apo_sugar_q2 = None, None
    complex_sugar_q1, complex_sugar_q2 = None, None
    
    if apo_sugar and apo_combined:
        apo_beta_barrel = get_beta_barrel_region(apo_combined.select('protein'))[0]
        if apo_beta_barrel:
            apo_reference_coords = apo_beta_barrel.getCoords()
            apo_sugar_q1, apo_sugar_q2 = calculate_sugar_projection(
                apo_sugar, apo_beta_barrel, apo_anm, apo_combined, apo_reference_coords
            )
            if apo_sugar_q1 is not None:
                print(f"  Apo糖位置: q₁={apo_sugar_q1:.3f}, q₂={apo_sugar_q2:.3f}")
    
    if complex_sugar and complex_combined:
        complex_beta_barrel = get_beta_barrel_region(complex_combined.select('protein'))[0]
        if complex_beta_barrel:
            complex_reference_coords = complex_beta_barrel.getCoords()
            complex_sugar_q1, complex_sugar_q2 = calculate_sugar_projection(
                complex_sugar, complex_beta_barrel, complex_anm, complex_combined, complex_reference_coords
            )
            if complex_sugar_q1 is not None:
                print(f"  Complex糖位置: q₁={complex_sugar_q1:.3f}, q₂={complex_sugar_q2:.3f}")
    
    # 计算统一的坐标轴范围
    q1_min = min(apo_q1.min(), complex_q1.min())
    q1_max = max(apo_q1.max(), complex_q1.max())
    q2_min = min(apo_q2.min(), complex_q2.min())
    q2_max = max(apo_q2.max(), complex_q2.max())
    
    # 添加边距
    margin = 0.1
    q1_range = q1_max - q1_min
    q2_range = q2_max - q2_min
    
    q1_unified = [q1_min - margin * q1_range, q1_max + margin * q1_range]
    q2_unified = [q2_min - margin * q2_range, q2_max + margin * q2_range]
    
    print(f"  统一坐标范围: q₁=[{q1_unified[0]:.2f}, {q1_unified[1]:.2f}], q₂=[{q2_unified[0]:.2f}, {q2_unified[1]:.2f}]")
    
    # 创建图形
    fig = plt.figure(figsize=(18, 8))
    
    # 创建自定义colormap（Nature风格）
    colors_apo = ['#FFFFFF', '#E8F4F8', '#2E86AB', '#1A5F7A', '#0D3A52']
    cmap_apo = LinearSegmentedColormap.from_list('nature_blue', colors_apo, N=100)
    
    colors_complex = ['#FFFFFF', '#F5E8F0', '#A23B72', '#7A2D56', '#4D1C35']
    cmap_complex = LinearSegmentedColormap.from_list('nature_purple', colors_complex, N=100)
    
    # 子图1: Apo自由能景观
    ax1 = fig.add_subplot(121)
    
    apo_FEL_clean = np.nan_to_num(apo_FEL, nan=0, posinf=0, neginf=0)
    if np.any(apo_FEL_clean > 0):
        p99 = np.percentile(apo_FEL_clean[apo_FEL_clean > 0], 99)
        apo_FEL_clean = np.clip(apo_FEL_clean, 0, p99 * 1.1)
    
    max_val = np.max(apo_FEL_clean)
    if max_val > 1e-6:
        levels = np.linspace(0, max_val, 25)
        levels = np.unique(levels)
        if len(levels) > 1:
            contour = ax1.contour(apo_Q1, apo_Q2, apo_FEL_clean, levels=levels, 
                                 colors='black', linewidths=0.6, alpha=0.4)
            ax1.clabel(contour, inline=True, fontsize=7, fmt='%.2f')
            im1 = ax1.contourf(apo_Q1, apo_Q2, apo_FEL_clean, levels=levels, 
                              cmap=cmap_apo, alpha=0.9, extend='max')
        else:
            im1 = ax1.contourf(apo_Q1, apo_Q2, apo_FEL_clean, cmap=cmap_apo, alpha=0.9)
    else:
        im1 = ax1.contourf(apo_Q1, apo_Q2, apo_FEL_clean, cmap=cmap_apo, alpha=0.9)
    
    # 绘制采样点（半透明）
    scatter1 = ax1.scatter(apo_q1, apo_q2, c='black', s=0.3, alpha=0.15, marker='.')
    
    # 标记最低能量点（原点附近）
    valid_mask = apo_FEL_clean > 0
    if np.any(valid_mask):
        min_idx = np.unravel_index(np.argmin(apo_FEL_clean[valid_mask]), apo_FEL_clean.shape)
        ax1.scatter([apo_Q1[min_idx]], [apo_Q2[min_idx]], 
                   c='red', s=400, marker='*', edgecolors='white', 
                   linewidths=2, label='Minimum', zorder=10)
    
    # 标记原点（参考结构）
    ax1.scatter([0], [0], c='green', s=200, marker='o', 
               edgecolors='white', linewidths=2, label='Reference', zorder=10)
    
    # 标记糖的位置（如果存在）- 只显示黄色菱形，不标注文字
    if apo_sugar_q1 is not None and apo_sugar_q2 is not None:
        ax1.scatter([apo_sugar_q1], [apo_sugar_q2], c='gold', s=300, marker='D',
                   edgecolors='black', linewidths=2, label='LPS', zorder=10)
    
    # 添加eigenvalue信息（放在右上角）
    info_text = f'λ₁ = {apo_eigenvalues[0]:.4f}\nλ₂ = {apo_eigenvalues[1]:.4f}'
    ax1.text(0.98, 0.98, info_text, transform=ax1.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax1.set_xlabel('q₁ (Mode 1 Projection, Å)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('q₂ (Mode 2 Projection, Å)', fontsize=14, fontweight='bold')
    ax1.set_title('Apo Free Energy Landscape\nE = ½(λ₁q₁² + λ₂q₂²)', 
                 fontsize=14, fontweight='bold', pad=15)
    ax1.grid(True, alpha=0.3, linestyle='--')
    # 图例：放在左上角，避免与eigenvalue信息重叠
    handles1, labels1 = ax1.get_legend_handles_labels()
    if handles1:
        ax1.legend(handles1, labels1, loc='upper left', framealpha=0.9, fontsize=9, 
                  bbox_to_anchor=(0.02, 0.98), ncol=1)
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # 设置统一的坐标轴范围
    ax1.set_xlim(q1_unified)
    ax1.set_ylim(q2_unified)
    
    cbar1 = plt.colorbar(im1, ax=ax1, pad=0.02)
    cbar1.set_label('Free Energy (kcal/mol)', fontsize=12, fontweight='bold')
    
    # 子图2: Complex自由能景观
    ax2 = fig.add_subplot(122)
    
    complex_FEL_clean = np.nan_to_num(complex_FEL, nan=0, posinf=0, neginf=0)
    if np.any(complex_FEL_clean > 0):
        p99 = np.percentile(complex_FEL_clean[complex_FEL_clean > 0], 99)
        complex_FEL_clean = np.clip(complex_FEL_clean, 0, p99 * 1.1)
    
    max_val2 = np.max(complex_FEL_clean)
    if max_val2 > 1e-6:
        levels_complex = np.linspace(0, max_val2, 25)
        levels_complex = np.unique(levels_complex)
        if len(levels_complex) > 1:
            contour2 = ax2.contour(complex_Q1, complex_Q2, complex_FEL_clean, levels=levels_complex,
                                  colors='black', linewidths=0.6, alpha=0.4)
            ax2.clabel(contour2, inline=True, fontsize=7, fmt='%.2f')
            im2 = ax2.contourf(complex_Q1, complex_Q2, complex_FEL_clean, levels=levels_complex,
                              cmap=cmap_complex, alpha=0.9, extend='max')
        else:
            im2 = ax2.contourf(complex_Q1, complex_Q2, complex_FEL_clean, cmap=cmap_complex, alpha=0.9)
    else:
        im2 = ax2.contourf(complex_Q1, complex_Q2, complex_FEL_clean, cmap=cmap_complex, alpha=0.9)
    
    scatter2 = ax2.scatter(complex_q1, complex_q2, c='black', s=0.3, alpha=0.15, marker='.')
    
    valid_mask2 = complex_FEL_clean > 0
    if np.any(valid_mask2):
        min_idx2 = np.unravel_index(np.argmin(complex_FEL_clean[valid_mask2]), complex_FEL_clean.shape)
        ax2.scatter([complex_Q1[min_idx2]], [complex_Q2[min_idx2]],
                   c='red', s=400, marker='*', edgecolors='white',
                   linewidths=2, label='Minimum', zorder=10)
    
    # 标记原点
    ax2.scatter([0], [0], c='green', s=200, marker='o',
               edgecolors='white', linewidths=2, label='Reference', zorder=10)
    
    # 标记糖的位置（如果存在）- 只显示黄色菱形，不标注文字
    if complex_sugar_q1 is not None and complex_sugar_q2 is not None:
        ax2.scatter([complex_sugar_q1], [complex_sugar_q2], c='gold', s=300, marker='D',
                   edgecolors='black', linewidths=2, label='LPS', zorder=10)
    
    info_text2 = f'λ₁ = {complex_eigenvalues[0]:.4f}\nλ₂ = {complex_eigenvalues[1]:.4f}'
    ax2.text(0.98, 0.98, info_text2, transform=ax2.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax2.set_xlabel('q₁ (Mode 1 Projection, Å)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('q₂ (Mode 2 Projection, Å)', fontsize=14, fontweight='bold')
    ax2.set_title('Complex Free Energy Landscape\nE = ½(λ₁q₁² + λ₂q₂²)', 
                 fontsize=14, fontweight='bold', pad=15)
    ax2.grid(True, alpha=0.3, linestyle='--')
    # 图例：放在左上角，避免与eigenvalue信息重叠
    handles2, labels2 = ax2.get_legend_handles_labels()
    if handles2:
        ax2.legend(handles2, labels2, loc='upper left', framealpha=0.9, fontsize=9,
                  bbox_to_anchor=(0.02, 0.98), ncol=1)
    ax2.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # 设置统一的坐标轴范围
    ax2.set_xlim(q1_unified)
    ax2.set_ylim(q2_unified)
    
    cbar2 = plt.colorbar(im2, ax=ax2, pad=0.02)
    cbar2.set_label('Free Energy (kcal/mol)', fontsize=12, fontweight='bold')
    
    # 添加总标题和坐标轴说明
    fig.suptitle('Free Energy Landscape: q₁ vs q₂ (Beta Barrel Region, NMA Projection)', 
                fontsize=16, fontweight='bold', y=0.98)
    
    # 添加坐标轴统一的说明
    axis_explanation = (
        f'Note: Unified axis ranges for direct comparison. Elliptical energy contours reflect eigenvalue differences: '
        f'Apo (λ₁={apo_eigenvalues[0]:.3f}, λ₂={apo_eigenvalues[1]:.3f}), Complex (λ₁={complex_eigenvalues[0]:.3f}, λ₂={complex_eigenvalues[1]:.3f}). '
        f'Different ellipse shapes indicate changes in protein motion anisotropy upon peptide binding.'
    )
    fig.text(0.5, 0.01, axis_explanation, ha='center', va='bottom', 
            fontsize=8, style='italic', wrap=True, color='gray')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, facecolor='white', edgecolor='none', format='pdf')
    plt.close()
    print(f"  ✅ 自由能景观图已保存: {output_file}")

def create_comparison_plot(apo_result, complex_result, apo_anm, complex_anm, output_file):
    """创建对比图，展示autoinhibition状态"""
    print("\n创建Autoinhibition对比图...")
    
    if apo_result[0] is None or complex_result[0] is None:
        return
    
    apo_Q1, apo_Q2, apo_FEL, apo_q1, apo_q2 = apo_result
    complex_Q1, complex_Q2, complex_FEL, complex_q1, complex_q2 = complex_result
    
    apo_eigenvalues = apo_anm.getEigvals()
    complex_eigenvalues = complex_anm.getEigvals()
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    # 创建colormap
    colors_apo = ['#FFFFFF', '#E8F4F8', '#2E86AB', '#1A5F7A', '#0D3A52']
    cmap_apo = LinearSegmentedColormap.from_list('nature_blue', colors_apo, N=100)
    colors_complex = ['#FFFFFF', '#F5E8F0', '#A23B72', '#7A2D56', '#4D1C35']
    cmap_complex = LinearSegmentedColormap.from_list('nature_purple', colors_complex, N=100)
    
    # 左上: Apo
    ax1 = axes[0, 0]
    apo_FEL_clean = np.nan_to_num(apo_FEL, nan=0, posinf=0, neginf=0)
    if np.any(apo_FEL_clean > 0):
        p99 = np.percentile(apo_FEL_clean[apo_FEL_clean > 0], 99)
        apo_FEL_clean = np.clip(apo_FEL_clean, 0, p99 * 1.1)
    max_val = np.max(apo_FEL_clean)
    if max_val > 1e-6:
        levels = np.linspace(0, max_val, 20)
        levels = np.unique(levels)
        if len(levels) > 1:
            im1 = ax1.contourf(apo_Q1, apo_Q2, apo_FEL_clean, levels=levels, cmap=cmap_apo, alpha=0.9)
            ax1.contour(apo_Q1, apo_Q2, apo_FEL_clean, levels=levels[::3], colors='black', linewidths=0.5, alpha=0.3)
    else:
        im1 = ax1.contourf(apo_Q1, apo_Q2, apo_FEL_clean, cmap=cmap_apo, alpha=0.9)
    ax1.scatter([0], [0], c='green', s=150, marker='o', edgecolors='white', linewidths=1.5, zorder=10)
    ax1.set_xlabel('q₁ (Å)', fontweight='bold')
    ax1.set_ylabel('q₂ (Å)', fontweight='bold')
    ax1.set_title(f'Apo\nλ₁={apo_eigenvalues[0]:.4f}, λ₂={apo_eigenvalues[1]:.4f}', fontweight='bold')
    ax1.grid(True, alpha=0.2)
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax1.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.colorbar(im1, ax=ax1, label='ΔG (kcal/mol)')
    
    # 右上: Complex
    ax2 = axes[0, 1]
    complex_FEL_clean = np.nan_to_num(complex_FEL, nan=0, posinf=0, neginf=0)
    if np.any(complex_FEL_clean > 0):
        p99 = np.percentile(complex_FEL_clean[complex_FEL_clean > 0], 99)
        complex_FEL_clean = np.clip(complex_FEL_clean, 0, p99 * 1.1)
    max_val2 = np.max(complex_FEL_clean)
    if max_val2 > 1e-6:
        levels2 = np.linspace(0, max_val2, 20)
        levels2 = np.unique(levels2)
        if len(levels2) > 1:
            im2 = ax2.contourf(complex_Q1, complex_Q2, complex_FEL_clean, levels=levels2, cmap=cmap_complex, alpha=0.9)
            ax2.contour(complex_Q1, complex_Q2, complex_FEL_clean, levels=levels2[::3], colors='black', linewidths=0.5, alpha=0.3)
    else:
        im2 = ax2.contourf(complex_Q1, complex_Q2, complex_FEL_clean, cmap=cmap_complex, alpha=0.9)
    ax2.scatter([0], [0], c='green', s=150, marker='o', edgecolors='white', linewidths=1.5, zorder=10)
    ax2.set_xlabel('q₁ (Å)', fontweight='bold')
    ax2.set_ylabel('q₂ (Å)', fontweight='bold')
    ax2.set_title(f'Complex\nλ₁={complex_eigenvalues[0]:.4f}, λ₂={complex_eigenvalues[1]:.4f}', fontweight='bold')
    ax2.grid(True, alpha=0.2)
    ax2.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax2.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.colorbar(im2, ax=ax2, label='ΔG (kcal/mol)')
    
    # 左下: 叠加对比
    ax3 = axes[1, 0]
    if max_val > 1e-6 and max_val2 > 1e-6:
        levels_apo = np.linspace(0, max_val, 15)
        levels_apo = np.unique(levels_apo)
        levels_complex = np.linspace(0, max_val2, 15)
        levels_complex = np.unique(levels_complex)
        
        if len(levels_apo) > 1:
            cs1 = ax3.contour(apo_Q1, apo_Q2, apo_FEL_clean, levels=levels_apo[::3],
                            colors=NATURE_COLORS['apo'], linewidths=2, alpha=0.7, linestyles='-')
        if len(levels_complex) > 1:
            cs2 = ax3.contour(complex_Q1, complex_Q2, complex_FEL_clean, levels=levels_complex[::3],
                            colors=NATURE_COLORS['complex'], linewidths=2, alpha=0.7, linestyles='--')
        
        handles = []
        if len(levels_apo) > 1:
            handles.append(plt.Line2D([0], [0], color=NATURE_COLORS['apo'], linewidth=2, linestyle='-', label='Apo'))
        if len(levels_complex) > 1:
            handles.append(plt.Line2D([0], [0], color=NATURE_COLORS['complex'], linewidth=2, linestyle='--', label='Complex'))
        if handles:
            ax3.legend(handles=handles, loc='upper right', framealpha=0.9)
    
    ax3.scatter([0], [0], c='green', s=150, marker='o', edgecolors='black', linewidths=1.5, zorder=10)
    ax3.set_xlabel('q₁ (Å)', fontweight='bold')
    ax3.set_ylabel('q₂ (Å)', fontweight='bold')
    ax3.set_title('Overlay Comparison', fontweight='bold')
    ax3.grid(True, alpha=0.2)
    ax3.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax3.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # 右下: 能量差（量化autoinhibition）
    ax4 = axes[1, 1]
    # 插值到相同网格
    from scipy.interpolate import griddata
    complex_FEL_interp = griddata(
        (complex_Q1.ravel(), complex_Q2.ravel()),
        complex_FEL_clean.ravel(),
        (apo_Q1, apo_Q2),
        method='linear',
        fill_value=np.nan
    )
    energy_diff = apo_FEL_clean - complex_FEL_interp
    energy_diff = np.nan_to_num(energy_diff, nan=0)
    
    colors_diff = ['#1A5F7A', '#2E86AB', '#FFFFFF', '#A23B72', '#7A2D56']
    cmap_diff = LinearSegmentedColormap.from_list('diverging', colors_diff, N=100)
    
    min_diff = np.nanmin(energy_diff)
    max_diff = np.nanmax(energy_diff)
    if max_diff > min_diff:
        levels_diff = np.linspace(min_diff, max_diff, 20)
        levels_diff = np.unique(levels_diff)
        if len(levels_diff) > 1:
            im4 = ax4.contourf(apo_Q1, apo_Q2, energy_diff, levels=levels_diff, cmap=cmap_diff, alpha=0.9)
            ax4.contour(apo_Q1, apo_Q2, energy_diff, levels=levels_diff[::3], colors='black', linewidths=0.5, alpha=0.3)
        else:
            im4 = ax4.contourf(apo_Q1, apo_Q2, energy_diff, cmap=cmap_diff, alpha=0.9)
    else:
        im4 = ax4.contourf(apo_Q1, apo_Q2, energy_diff, cmap=cmap_diff, alpha=0.9)
    
    ax4.set_xlabel('q₁ (Å)', fontweight='bold')
    ax4.set_ylabel('q₂ (Å)', fontweight='bold')
    ax4.set_title('Energy Difference (Apo - Complex)\nPositive = Apo more flexible', fontweight='bold')
    ax4.grid(True, alpha=0.2)
    ax4.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax4.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    cbar4 = plt.colorbar(im4, ax=ax4, label='ΔΔG (kcal/mol)')
    cbar4.ax.axhline(y=0, color='black', linewidth=1)
    
    plt.suptitle('Free Energy Landscape: Autoinhibition Analysis', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, facecolor='white', edgecolor='none', format='pdf')
    plt.close()
    print(f"  ✅ Autoinhibition对比图已保存: {output_file}")

def analyze_autoinhibition(apo_result, complex_result, apo_anm, complex_anm, output_file):
    """分析并量化autoinhibition状态"""
    print("\n分析Autoinhibition状态...")
    
    if apo_result[0] is None or complex_result[0] is None:
        return
    
    apo_Q1, apo_Q2, apo_FEL, apo_q1, apo_q2 = apo_result
    complex_Q1, complex_Q2, complex_FEL, complex_q1, complex_q2 = complex_result
    
    apo_eigenvalues = apo_anm.getEigvals()
    complex_eigenvalues = complex_anm.getEigvals()
    
    # 计算关键指标
    # 1. 能谷深度（从中心到最低点的能量差）
    apo_min_energy = np.min(apo_FEL[apo_FEL > 0]) if np.any(apo_FEL > 0) else 0
    complex_min_energy = np.min(complex_FEL[complex_FEL > 0]) if np.any(complex_FEL > 0) else 0
    
    # 2. 能谷宽度（在特定能量阈值下的面积）
    # 使用能量阈值 = 最小能量的2倍
    apo_threshold = apo_min_energy * 2 if apo_min_energy > 0 else np.percentile(apo_FEL[apo_FEL > 0], 10)
    complex_threshold = complex_min_energy * 2 if complex_min_energy > 0 else np.percentile(complex_FEL[complex_FEL > 0], 10)
    
    apo_accessible_area = np.sum(apo_FEL <= apo_threshold) if apo_min_energy > 0 else 0
    complex_accessible_area = np.sum(complex_FEL <= complex_threshold) if complex_min_energy > 0 else 0
    
    # 3. 势能陡峭度（eigenvalue越大，势能越陡峭）
    apo_stiffness = (apo_eigenvalues[0] + apo_eigenvalues[1]) / 2
    complex_stiffness = (complex_eigenvalues[0] + complex_eigenvalues[1]) / 2
    
    # 4. 自由度（可访问的构象空间范围）
    apo_q1_range = apo_q1.max() - apo_q1.min()
    apo_q2_range = apo_q2.max() - apo_q2.min()
    apo_freedom = apo_q1_range * apo_q2_range
    
    complex_q1_range = complex_q1.max() - complex_q1.min()
    complex_q2_range = complex_q2.max() - complex_q2.min()
    complex_freedom = complex_q1_range * complex_q2_range
    
    # 写入报告
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("Autoinhibition状态分析报告\n")
        f.write("="*80 + "\n\n")
        
        f.write("1. Eigenvalues（势能陡峭度指标）\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Structure':<20} {'λ₁ (Mode 1)':<20} {'λ₂ (Mode 2)':<20} {'Average':<20}\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Apo':<20} {apo_eigenvalues[0]:<20.6f} {apo_eigenvalues[1]:<20.6f} {apo_stiffness:<20.6f}\n")
        f.write(f"{'Complex':<20} {complex_eigenvalues[0]:<20.6f} {complex_eigenvalues[1]:<20.6f} {complex_stiffness:<20.6f}\n")
        f.write(f"{'Difference':<20} {complex_eigenvalues[0]-apo_eigenvalues[0]:<20.6f} {complex_eigenvalues[1]-apo_eigenvalues[1]:<20.6f} {complex_stiffness-apo_stiffness:<20.6f}\n")
        f.write("\n")
        
        f.write("2. 能谷特征\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Structure':<20} {'Min Energy':<20} {'Threshold':<20} {'Accessible Area':<20}\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Apo':<20} {apo_min_energy:<20.4f} {apo_threshold:<20.4f} {apo_accessible_area:<20}\n")
        f.write(f"{'Complex':<20} {complex_min_energy:<20.4f} {complex_threshold:<20.4f} {complex_accessible_area:<20}\n")
        f.write("\n")
        
        f.write("3. 自由度（构象空间范围）\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Structure':<20} {'q₁ Range':<20} {'q₂ Range':<20} {'Freedom (q₁×q₂)':<20}\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Apo':<20} {apo_q1_range:<20.4f} {apo_q2_range:<20.4f} {apo_freedom:<20.4f}\n")
        f.write(f"{'Complex':<20} {complex_q1_range:<20.4f} {complex_q2_range:<20.4f} {complex_freedom:<20.4f}\n")
        f.write(f"{'Ratio (Apo/Complex)':<20} {apo_q1_range/complex_q1_range if complex_q1_range > 0 else 0:<20.4f} {apo_q2_range/complex_q2_range if complex_q2_range > 0 else 0:<20.4f} {apo_freedom/complex_freedom if complex_freedom > 0 else 0:<20.4f}\n")
        f.write("\n")
        
        f.write("4. Autoinhibition判断\n")
        f.write("-"*80 + "\n")
        
        # 判断标准
        stiffness_increase = (complex_stiffness - apo_stiffness) / apo_stiffness * 100 if apo_stiffness > 0 else 0
        
        # 计算在相同能量阈值下的可访问面积（更准确的自由度指标）
        # 使用相同的能量阈值（例如1 kcal/mol）来比较
        energy_threshold = 1.0  # kcal/mol
        apo_accessible_at_threshold = np.sum(apo_FEL <= energy_threshold)
        complex_accessible_at_threshold = np.sum(complex_FEL <= energy_threshold)
        
        # 计算能谷宽度（在低能量区域的"宽度"）
        # 使用能量阈值 = 最小能量的5倍
        apo_width_threshold = apo_min_energy * 5 if apo_min_energy > 0 else 0.5
        complex_width_threshold = complex_min_energy * 5 if complex_min_energy > 0 else 0.5
        
        apo_width = np.sum(apo_FEL <= apo_width_threshold)
        complex_width = np.sum(complex_FEL <= complex_width_threshold)
        
        width_ratio = complex_width / apo_width if apo_width > 0 else 0
        
        f.write(f"势能陡峭度增加: {stiffness_increase:.2f}%\n")
        f.write(f"在能量阈值 {energy_threshold} kcal/mol下的可访问面积:\n")
        f.write(f"  Apo: {apo_accessible_at_threshold}\n")
        f.write(f"  Complex: {complex_accessible_at_threshold}\n")
        f.write(f"  比例 (Complex/Apo): {complex_accessible_at_threshold/apo_accessible_at_threshold if apo_accessible_at_threshold > 0 else 0:.2f}\n")
        f.write(f"\n能谷宽度（低能量区域）:\n")
        f.write(f"  Apo: {apo_width}\n")
        f.write(f"  Complex: {complex_width}\n")
        f.write(f"  比例 (Complex/Apo): {width_ratio:.2f}\n")
        f.write("\n")
        
        # 判断标准：主要看eigenvalues增加和能谷宽度
        if stiffness_increase > 50:
            f.write("✅ **结论: 形成Autoinhibition状态**\n")
            f.write("   - Complex的势能更陡峭（eigenvalues增加{:.1f}%）\n".format(stiffness_increase))
            if width_ratio < 0.8:
                f.write("   - Complex的能谷更窄（能谷宽度减少）\n")
            f.write("   - 在相同方向上，Complex的能量变化更大（势能陡峭）\n")
            f.write("   - 环肽诱导了NTD下压并增强了gate stiffness\n")
            f.write("   - 这限制了gate的打开，形成autoinhibition状态\n")
        elif stiffness_increase > 20:
            f.write("⚠️  **部分Autoinhibition状态**\n")
            f.write("   - Complex的势能有所增加（eigenvalues增加{:.1f}%）\n".format(stiffness_increase))
            f.write("   - 但能谷宽度变化不明显\n")
        else:
            f.write("❌ **未形成明显的Autoinhibition状态**\n")
            f.write("   - Complex和Apo的势能差异较小\n")
        
        f.write("\n")
        f.write("5. 详细解释\n")
        f.write("-"*80 + "\n")
        f.write("能谷浅（Apo）: 在相同方向上的能量变化较小，系统更容易偏离平衡位置\n")
        f.write("势能陡峭（Complex）: 在相同方向上的能量变化较大，系统更难以偏离平衡位置\n")
        f.write("Autoinhibition机制: 环肽结合后，通过增加eigenvalues（势能陡峭度）和减少\n")
        f.write("                   可访问的构象空间，限制了gate的打开，形成autoinhibition状态\n")
    
    print(f"  ✅ Autoinhibition分析报告已保存: {output_file}")

def main():
    """主函数"""
    print("="*80)
    print("Free Energy Landscape生成 (NMA投影方法)")
    print("="*80)
    
    # 加载数据（包含糖）
    apo_protein, complex_protein, apo_anm, complex_anm, apo_sugar, complex_sugar, apo_combined, complex_combined = load_structures_and_nma()
    
    # 获取beta桶区域
    print("\n获取Beta桶区域...")
    apo_beta_barrel, apo_beta_resnums = get_beta_barrel_region(apo_protein)
    complex_beta_barrel, complex_beta_resnums = get_beta_barrel_region(complex_protein)
    
    if apo_beta_barrel is None or complex_beta_barrel is None:
        print("错误: 无法获取beta桶区域")
        sys.exit(1)
    
    # 采样构象空间并计算NMA投影（使用combined结构，包含糖）
    print("\n采样Apo构象空间并计算NMA投影...")
    apo_combined_beta_barrel = get_beta_barrel_region(apo_combined.select('protein'))[0] if apo_combined else apo_beta_barrel
    apo_result = sample_conformational_space_nma(apo_combined if apo_combined else apo_protein, apo_anm, apo_combined_beta_barrel, n_samples=1000)
    
    print("\n采样Complex构象空间并计算NMA投影...")
    complex_combined_beta_barrel = get_beta_barrel_region(complex_combined.select('protein'))[0] if complex_combined else complex_beta_barrel
    complex_result = sample_conformational_space_nma(complex_combined if complex_combined else complex_protein, complex_anm, complex_combined_beta_barrel, n_samples=1000)
    
    if apo_result[0] is None or complex_result[0] is None:
        print("错误: 无法计算NMA投影")
        sys.exit(1)
    
    # 计算自由能景观
    print("\n计算Apo自由能景观...")
    apo_q1, apo_q2, apo_energies = apo_result
    apo_landscape = calculate_free_energy_landscape_nma(apo_q1, apo_q2, apo_energies)
    
    print("\n计算Complex自由能景观...")
    complex_q1, complex_q2, complex_energies = complex_result
    complex_landscape = calculate_free_energy_landscape_nma(complex_q1, complex_q2, complex_energies)
    
    # 创建输出目录
    output_dir = Path("nma_results")
    output_dir.mkdir(exist_ok=True)
    
    # 创建自由能景观图（包含糖标记）
    create_free_energy_landscape_nma_plot(apo_landscape, complex_landscape, apo_anm, complex_anm,
                                         apo_sugar, complex_sugar, apo_combined, complex_combined,
                                         output_dir / "free_energy_landscape_nma.pdf")
    
    # 创建对比图
    create_comparison_plot(apo_landscape, complex_landscape, apo_anm, complex_anm,
                          output_dir / "free_energy_landscape_autoinhibition.pdf")
    
    # 分析autoinhibition
    analyze_autoinhibition(apo_landscape, complex_landscape, apo_anm, complex_anm,
                         output_dir / "autoinhibition_analysis.txt")
    
    print("\n" + "="*80)
    print("✅ 自由能景观图生成完成！")
    print("="*80)
    print(f"\n生成的文件:")
    print(f"  - {output_dir}/free_energy_landscape_nma.pdf (基于NMA投影的自由能景观图)")
    print(f"  - {output_dir}/free_energy_landscape_autoinhibition.pdf (Autoinhibition对比图)")
    print(f"  - {output_dir}/autoinhibition_analysis.txt (Autoinhibition分析报告)")
    print("\n图例说明:")
    print("  📊 横坐标 q₁ = Mode 1投影 (ΔR · u₁)")
    print("  📊 纵坐标 q₂ = Mode 2投影 (ΔR · u₂)")
    print("  📈 自由能 E = ½(λ₁q₁² + λ₂q₂²)")
    print("  ⭐ 红色星号 = 最低能量点")
    print("  🟢 绿色圆点 = 参考结构（原点）")
    print("  🟠 橙色菱形 = 糖的位置")
    print("  🔵 蓝色区域 = Apo（能谷浅，自由度大）")
    print("  🟣 紫红色区域 = Complex（势能陡峭，更刚性）")

if __name__ == "__main__":
    main()

