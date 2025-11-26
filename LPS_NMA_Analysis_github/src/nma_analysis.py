#!/usr/bin/env python3
"""
ProDy NMA分析脚本
比较apo和complex结构，分析环肽对NTD和桶蛋白打开的影响
"""

import prody
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

# 设置matplotlib为Nature风格
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'axes.linewidth': 1.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1
})

def load_structure(pdb_file, name):
    """加载PDB结构，只保留蛋白质部分"""
    print(f"\n{'='*60}")
    print(f"加载结构: {name}")
    print(f"{'='*60}")
    
    # 解析PDB文件
    structure = prody.parsePDB(pdb_file)
    
    # 选择蛋白质（排除配体和水）
    protein = structure.select('protein')
    
    if protein is None:
        raise ValueError(f"无法从 {pdb_file} 中提取蛋白质")
    
    print(f"  总原子数: {structure.numAtoms()}")
    print(f"  蛋白质原子数: {protein.numAtoms()}")
    print(f"  残基数: {len(set(protein.getResnums()))}")
    print(f"  链数: {len(set(protein.getChids()))}")
    
    return protein

def perform_nma(protein, name, n_modes=20):
    """执行正常模式分析"""
    print(f"\n{'='*60}")
    print(f"执行NMA: {name}")
    print(f"{'='*60}")
    
    # 使用ANM (Anisotropic Network Model)
    print("  使用ANM模型...")
    anm = prody.ANM(f'{name}_ANM')
    
    # 构建Hessian矩阵
    print("  构建Hessian矩阵...")
    anm.buildHessian(protein)
    
    # 计算正常模式
    print(f"  计算前{n_modes}个正常模式...")
    anm.calcModes(n_modes=n_modes)
    
    # 获取eigenvalues
    eigenvalues = anm.getEigvals()
    
    print(f"  完成！前6个eigenvalues:")
    for i in range(min(6, len(eigenvalues))):
        print(f"    Mode {i+1}: {eigenvalues[i]:.6f}")
    
    return anm, eigenvalues

def calculate_rms_displacement(anm, mode_indices=[0, 1, 2, 3, 4, 5]):
    """计算RMS displacement"""
    rms_displacements = []
    
    for mode_idx in mode_indices:
        if mode_idx >= len(anm):
            break
        
        # 获取模式向量
        mode = anm[mode_idx]
        
        # 计算每个原子的位移幅度
        displacements = mode.getArrayNx3()
        rms_per_atom = np.sqrt(np.sum(displacements**2, axis=1))
        rms_mean = np.mean(rms_per_atom)
        
        rms_displacements.append(rms_mean)
    
    return rms_displacements

def identify_key_residues(protein, name=""):
    """识别关键残基：NTD、β1、β26"""
    key_residues = {
        'NTD': [],
        'beta1': [],
        'beta26': []
    }
    
    # 获取链信息
    chains = sorted(set(protein.getChids()))
    print(f"\n  检测到链: {chains}")
    
    # 分析每个链的长度
    chain_info = {}
    for chid in chains:
        chain = protein.select(f'chain {chid}')
        if chain:
            resnums = sorted(set(chain.getResnums()))
            chain_info[chid] = {
                'resnums': resnums,
                'start': min(resnums),
                'end': max(resnums),
                'length': len(resnums)
            }
            print(f"    链 {chid}: 残基 {min(resnums)} - {max(resnums)} ({len(resnums)} 个残基)")
    
    # NTD: 通常是最短的链或N端区域
    # 假设最短的链是NTD
    if chain_info:
        shortest_chain = min(chain_info.items(), key=lambda x: x[1]['length'])
        ntd_chain = shortest_chain[0]
        ntd_residues = chain_info[ntd_chain]['resnums']
        key_residues['NTD'] = ntd_residues
        print(f"  → NTD识别为链 {ntd_chain}: 残基 {min(ntd_residues)} - {max(ntd_residues)}")
    
    # β1和β26: 假设是主链（最长的链）中的特定区域
    # 这里需要根据实际桶蛋白结构来定义
    # 对于桶蛋白，β链通常是连续的残基段
    if chain_info:
        longest_chain = max(chain_info.items(), key=lambda x: x[1]['length'])
        main_chain = longest_chain[0]
        main_resnums = chain_info[main_chain]['resnums']
        
        # β1: 假设是主链的前1/4区域（需要根据实际结构调整）
        beta1_start = main_resnums[0]
        beta1_end = main_resnums[len(main_resnums)//4]
        beta1_residues = [r for r in main_resnums if beta1_start <= r <= beta1_end]
        key_residues['beta1'] = beta1_residues[:20]  # 取前20个作为β1
        print(f"  → β1识别为链 {main_chain} 的前段: 残基 {beta1_start} - {beta1_end} (取前20个)")
        
        # β26: 假设是主链的后1/4区域（需要根据实际结构调整）
        beta26_start = main_resnums[3*len(main_resnums)//4]
        beta26_end = main_resnums[-1]
        beta26_residues = [r for r in main_resnums if beta26_start <= r <= beta26_end]
        key_residues['beta26'] = beta26_residues[-20:]  # 取后20个作为β26
        print(f"  → β26识别为链 {main_chain} 的后段: 残基 {beta26_start} - {beta26_end} (取后20个)")
    
    print("\n  ⚠️  注意: 以上识别为自动推测，请根据实际结构验证！")
    print("    如需调整，请修改identify_key_residues函数或手动指定残基编号")
    
    return key_residues

def analyze_key_residue_motion(anm, protein, key_residues, name):
    """分析关键残基的运动幅度"""
    print(f"\n分析关键残基运动: {name}")
    
    # 获取残基编号
    resnums = protein.getResnums()
    
    motion_data = {}
    
    for region, residue_list in key_residues.items():
        if not residue_list:
            continue
        
        # 选择该区域的原子
        # ProDy的resnum选择语法需要逐个指定
        selection = None
        for resnum in residue_list:
            res_selection = protein.select(f'resnum {resnum}')
            if res_selection is not None:
                if selection is None:
                    selection = res_selection
                else:
                    selection = selection + res_selection
        
        if selection is None or selection.numAtoms() == 0:
            print(f"  警告: 无法选择{region}区域")
            continue
        
        # 获取这些原子的索引（相对于protein选择）
        # 需要找到selection中的原子在protein中的位置
        selection_indices = selection.getIndices()
        protein_indices = protein.getIndices()
        
        # 找到selection中的原子在protein中的对应索引
        indices_in_protein = []
        for sel_idx in selection_indices:
            if sel_idx in protein_indices:
                # 找到在protein中的位置
                idx_in_protein = np.where(protein_indices == sel_idx)[0]
                if len(idx_in_protein) > 0:
                    indices_in_protein.append(idx_in_protein[0])
        
        if not indices_in_protein:
            print(f"  警告: 无法找到{region}区域在protein中的索引")
            continue
        
        indices_in_protein = np.array(indices_in_protein)
        
        # 计算前6个模式的运动幅度
        motions = []
        for mode_idx in range(min(6, len(anm))):
            mode = anm[mode_idx]
            displacements = mode.getArrayNx3()
            
            # 确保索引在范围内
            valid_indices = indices_in_protein[indices_in_protein < len(displacements)]
            if len(valid_indices) == 0:
                continue
            
            # 计算该区域的平均RMS displacement
            region_displacements = displacements[valid_indices]
            rms = np.sqrt(np.sum(region_displacements**2, axis=1))
            avg_rms = np.mean(rms)
            motions.append(avg_rms)
        
        motion_data[region] = motions
        print(f"  {region}: 平均RMS displacement = {np.mean(motions):.4f} Å")
    
    return motion_data

def calculate_gate_openness(anm, protein, gate_residues_1, gate_residues_2):
    """
    计算gate openness（门开放度）
    通过测量两个gate区域之间的距离变化
    """
    print("\n计算Gate Openness...")
    
    # 获取gate区域的原子
    gate1_atoms = None
    for resnum in gate_residues_1:
        res_atoms = protein.select(f'resnum {resnum}')
        if res_atoms is not None:
            if gate1_atoms is None:
                gate1_atoms = res_atoms
            else:
                gate1_atoms = gate1_atoms + res_atoms
    
    gate2_atoms = None
    for resnum in gate_residues_2:
        res_atoms = protein.select(f'resnum {resnum}')
        if res_atoms is not None:
            if gate2_atoms is None:
                gate2_atoms = res_atoms
            else:
                gate2_atoms = gate2_atoms + res_atoms
    
    if gate1_atoms is None or gate2_atoms is None:
        print("  警告: 无法选择gate区域，跳过gate openness计算")
        return None
    
    # 计算初始距离（质心距离）
    gate1_center = gate1_atoms.getCoords().mean(axis=0)
    gate2_center = gate2_atoms.getCoords().mean(axis=0)
    initial_distance = np.linalg.norm(gate1_center - gate2_center)
    
    print(f"  初始gate距离: {initial_distance:.2f} Å")
    
    # 获取protein的索引
    protein_indices = protein.getIndices()
    gate1_selection_indices = gate1_atoms.getIndices()
    gate2_selection_indices = gate2_atoms.getIndices()
    
    # 计算前6个模式对gate距离的影响
    gate_openness = []
    
    for mode_idx in range(min(6, len(anm))):
        mode = anm[mode_idx]
        displacements = mode.getArrayNx3()
        
        # 找到gate区域在protein中的索引
        gate1_indices_in_protein = []
        for idx in gate1_selection_indices:
            if idx in protein_indices:
                pos = np.where(protein_indices == idx)[0]
                if len(pos) > 0:
                    gate1_indices_in_protein.append(pos[0])
        
        gate2_indices_in_protein = []
        for idx in gate2_selection_indices:
            if idx in protein_indices:
                pos = np.where(protein_indices == idx)[0]
                if len(pos) > 0:
                    gate2_indices_in_protein.append(pos[0])
        
        if not gate1_indices_in_protein or not gate2_indices_in_protein:
            continue
        
        gate1_indices_in_protein = np.array(gate1_indices_in_protein)
        gate2_indices_in_protein = np.array(gate2_indices_in_protein)
        
        # 确保索引在范围内
        gate1_valid = gate1_indices_in_protein[gate1_indices_in_protein < len(displacements)]
        gate2_valid = gate2_indices_in_protein[gate2_indices_in_protein < len(displacements)]
        
        if len(gate1_valid) == 0 or len(gate2_valid) == 0:
            continue
        
        gate1_disp = displacements[gate1_valid].mean(axis=0)
        gate2_disp = displacements[gate2_valid].mean(axis=0)
        
        # 计算位移方向对gate距离的影响
        # 如果两个gate区域沿连接方向移动，距离会增加（开放）
        gate_vector = gate2_center - gate1_center
        if np.linalg.norm(gate_vector) < 1e-6:
            continue
        gate_vector_norm = gate_vector / np.linalg.norm(gate_vector)
        
        # 投影位移到gate方向
        gate1_proj = np.dot(gate1_disp, -gate_vector_norm)  # 负方向（gate1向外）
        gate2_proj = np.dot(gate2_disp, gate_vector_norm)   # 正方向（gate2向外）
        
        # 总开放度 = 两个gate向外移动的总和
        openness = gate1_proj + gate2_proj
        gate_openness.append(openness)
    
    if not gate_openness:
        print("  警告: 无法计算gate openness")
        return None
    
    return {
        'initial_distance': initial_distance,
        'openness_per_mode': gate_openness,
        'average_openness': np.mean(gate_openness),
        'max_openness': np.max(gate_openness)
    }
def generate_nma_frames(anm, protein, name, n_frames=20, amplitude=5.0):
    """生成NMA动画帧（PDB文件）"""
    print(f"\n生成NMA动画帧: {name}")
    
    output_dir = Path(f"nma_frames_{name}")
    output_dir.mkdir(exist_ok=True)
    
    # 获取前6个低频模式
    n_modes_to_use = min(6, len(anm))
    
    for mode_idx in range(n_modes_to_use):
        mode = anm[mode_idx]
        
        # 生成多个振幅的帧
        for frame_idx, amp in enumerate(np.linspace(-amplitude, amplitude, n_frames)):
            # 计算位移
            displacements = mode.getArrayNx3() * amp
            
            # 应用位移
            new_coords = protein.getCoords() + displacements
            
            # 创建新的AtomGroup
            new_protein = protein.copy()
            new_protein.setCoords(new_coords)
            
            # 保存PDB
            filename = output_dir / f"mode{mode_idx+1}_frame{frame_idx:03d}_amp{amp:.2f}.pdb"
            prody.writePDB(str(filename), new_protein)
        
        print(f"  Mode {mode_idx+1}: 生成 {n_frames} 帧")
    
    print(f"  所有帧保存在: {output_dir}/")
    return output_dir

def plot_comparison(apo_anm, complex_anm, apo_eigenvalues, complex_eigenvalues, 
                   apo_motion, complex_motion, output_dir):
    """生成对比图"""
    print("\n生成对比图...")
    
    # 1. Eigenvalues对比
    fig, ax = plt.subplots(figsize=(8, 6))
    modes = np.arange(1, 7)
    ax.plot(modes, apo_eigenvalues[:6], 'o-', label='Apo', linewidth=2, markersize=8, color='#2E86AB')
    ax.plot(modes, complex_eigenvalues[:6], 's-', label='Complex', linewidth=2, markersize=8, color='#A23B72')
    ax.set_xlabel('Mode Number', fontsize=14, fontweight='bold')
    ax.set_ylabel('Eigenvalue', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12, frameon=False)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / 'eigenvalues_comparison.png', dpi=300)
    plt.close()
    print("  ✅ Eigenvalues对比图已保存")
    
    # 2. RMS Displacement对比
    apo_rms = calculate_rms_displacement(apo_anm)
    complex_rms = calculate_rms_displacement(complex_anm)
    
    fig, ax = plt.subplots(figsize=(8, 6))
    modes = np.arange(1, len(apo_rms)+1)
    width = 0.35
    x = np.arange(len(modes))
    ax.bar(x - width/2, apo_rms, width, label='Apo', color='#2E86AB', alpha=0.8)
    ax.bar(x + width/2, complex_rms, width, label='Complex', color='#A23B72', alpha=0.8)
    ax.set_xlabel('Mode Number', fontsize=14, fontweight='bold')
    ax.set_ylabel('RMS Displacement (Å)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(modes)
    ax.legend(fontsize=12, frameon=False)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig(output_dir / 'rms_displacement_comparison.png', dpi=300)
    plt.close()
    print("  ✅ RMS Displacement对比图已保存")
    
    # 3. 关键残基运动幅度对比
    if apo_motion and complex_motion:
        regions = list(set(apo_motion.keys()) | set(complex_motion.keys()))
        if regions:
            fig, axes = plt.subplots(len(regions), 1, figsize=(10, 4*len(regions)))
            if len(regions) == 1:
                axes = [axes]
            
            for idx, region in enumerate(regions):
                ax = axes[idx]
                apo_data = apo_motion.get(region, [])
                complex_data = complex_motion.get(region, [])
                
                if apo_data and complex_data:
                    modes = np.arange(1, min(len(apo_data), len(complex_data)) + 1)
                    width = 0.35
                    x = np.arange(len(modes))
                    ax.bar(x - width/2, apo_data[:len(modes)], width, label='Apo', 
                          color='#2E86AB', alpha=0.8)
                    ax.bar(x + width/2, complex_data[:len(modes)], width, label='Complex', 
                          color='#A23B72', alpha=0.8)
                    ax.set_xlabel('Mode Number', fontsize=12, fontweight='bold')
                    ax.set_ylabel('RMS Displacement (Å)', fontsize=12, fontweight='bold')
                    ax.set_title(f'{region} 运动幅度', fontsize=14, fontweight='bold')
                    ax.set_xticks(x)
                    ax.set_xticklabels(modes)
                    ax.legend(fontsize=10, frameon=False)
                    ax.grid(True, alpha=0.3, axis='y')
            
            plt.tight_layout()
            plt.savefig(output_dir / 'key_residue_motion_comparison.png', dpi=300)
            plt.close()
            print("  ✅ 关键残基运动幅度对比图已保存")

def generate_summary_table(apo_eigenvalues, complex_eigenvalues, apo_gate, complex_gate, output_file):
    """生成数值化表格"""
    print("\n生成数值化表格...")
    
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("NMA分析结果汇总表\n")
        f.write("="*80 + "\n\n")
        
        # Eigenvalues
        f.write("前6个Eigenvalues:\n")
        f.write("-"*80 + "\n")
        f.write(f"{'Mode':<10} {'Apo':<20} {'Complex':<20} {'Difference':<20}\n")
        f.write("-"*80 + "\n")
        for i in range(6):
            diff = complex_eigenvalues[i] - apo_eigenvalues[i]
            f.write(f"{i+1:<10} {apo_eigenvalues[i]:<20.6f} {complex_eigenvalues[i]:<20.6f} {diff:<20.6f}\n")
        f.write("\n")
        
        # Gate Openness
        if apo_gate and complex_gate:
            f.write("Gate Openness分析:\n")
            f.write("-"*80 + "\n")
            f.write(f"{'Metric':<30} {'Apo':<20} {'Complex':<20}\n")
            f.write("-"*80 + "\n")
            f.write(f"{'Average Openness':<30} {apo_gate['average_openness']:<20.4f} {complex_gate['average_openness']:<20.4f}\n")
            f.write(f"{'Max Openness':<30} {apo_gate['max_openness']:<20.4f} {complex_gate['max_openness']:<20.4f}\n")
            f.write(f"{'Initial Gate Distance (Å)':<30} {apo_gate['initial_distance']:<20.2f} {complex_gate['initial_distance']:<20.2f}\n")
            f.write("\n")
            
            # 趋势判断
            f.write("Autoinhibition趋势判断:\n")
            f.write("-"*80 + "\n")
            if apo_gate['average_openness'] > 0:
                f.write("Apo: 显示开放趋势 (Average Openness > 0)\n")
            else:
                f.write("Apo: 显示封闭趋势 (Average Openness < 0)\n")
            
            if complex_gate['average_openness'] < 0:
                f.write("Complex: 显示封闭趋势 (Average Openness < 0)\n")
            else:
                f.write("Complex: 显示开放趋势 (Average Openness > 0)\n")
            f.write("\n")
    
    print(f"  ✅ 表格已保存: {output_file}")

def main():
    """主函数"""
    print("="*80)
    print("ProDy NMA分析: Apo vs Complex")
    print("="*80)
    
    # 文件路径
    apo_pdb = "no-peptide.pdb"
    complex_pdb = "complex.pdb"
    
    # 检查文件
    if not os.path.exists(apo_pdb):
        print(f"错误: 找不到文件 {apo_pdb}")
        sys.exit(1)
    if not os.path.exists(complex_pdb):
        print(f"错误: 找不到文件 {complex_pdb}")
        sys.exit(1)
    
    # 创建输出目录
    output_dir = Path("nma_results")
    output_dir.mkdir(exist_ok=True)
    
    # 1. 加载结构
    apo_protein = load_structure(apo_pdb, "Apo")
    complex_protein = load_structure(complex_pdb, "Complex")
    
    # 2. 执行NMA
    apo_anm, apo_eigenvalues = perform_nma(apo_protein, "Apo", n_modes=20)
    complex_anm, complex_eigenvalues = perform_nma(complex_protein, "Complex", n_modes=20)
    
    # 3. 识别关键残基（需要用户根据实际结构调整）
    print("\n" + "="*60)
    print("关键残基识别")
    print("="*60)
    print("⚠️  请根据您的PDB结构手动定义以下残基:")
    print("   - NTD: N端结构域残基编号")
    print("   - β1: 第一条β链残基编号")
    print("   - β26: 第26条β链残基编号")
    print("   - Gate区域1和2: 用于计算gate openness的残基")
    
    # 识别关键残基
    apo_key_residues = identify_key_residues(apo_protein, "Apo")
    complex_key_residues = identify_key_residues(complex_protein, "Complex")
    
    # 4. 分析关键残基运动
    apo_motion = analyze_key_residue_motion(apo_anm, apo_protein, apo_key_residues, "Apo")
    complex_motion = analyze_key_residue_motion(complex_anm, complex_protein, complex_key_residues, "Complex")
    
    # 5. 计算Gate Openness（需要定义gate残基）
    print("\n" + "="*60)
    print("Gate Openness计算")
    print("="*60)
    print("⚠️  需要手动定义gate区域的残基编号")
    print("   请修改脚本中的gate_residues_1和gate_residues_2")
    
    # 自动识别gate区域
    # Gate通常是桶蛋白开口处的残基，可能是β1和β26附近的残基
    print("\n自动识别Gate区域...")
    apo_gate_residues_1 = apo_key_residues.get('beta1', [])[:10]  # β1的前10个残基
    apo_gate_residues_2 = apo_key_residues.get('beta26', [])[-10:]  # β26的后10个残基
    complex_gate_residues_1 = complex_key_residues.get('beta1', [])[:10]
    complex_gate_residues_2 = complex_key_residues.get('beta26', [])[-10:]
    
    if apo_gate_residues_1 and apo_gate_residues_2:
        print(f"  Apo Gate区域1: 残基 {min(apo_gate_residues_1)} - {max(apo_gate_residues_1)}")
        print(f"  Apo Gate区域2: 残基 {min(apo_gate_residues_2)} - {max(apo_gate_residues_2)}")
    if complex_gate_residues_1 and complex_gate_residues_2:
        print(f"  Complex Gate区域1: 残基 {min(complex_gate_residues_1)} - {max(complex_gate_residues_1)}")
        print(f"  Complex Gate区域2: 残基 {min(complex_gate_residues_2)} - {max(complex_gate_residues_2)}")
    
    apo_gate = None
    complex_gate = None
    
    if apo_gate_residues_1 and apo_gate_residues_2:
        apo_gate = calculate_gate_openness(apo_anm, apo_protein, apo_gate_residues_1, apo_gate_residues_2)
    if complex_gate_residues_1 and complex_gate_residues_2:
        complex_gate = calculate_gate_openness(complex_anm, complex_protein, complex_gate_residues_1, complex_gate_residues_2)
    
    # 6. 生成对比图
    plot_comparison(apo_anm, complex_anm, apo_eigenvalues, complex_eigenvalues, 
                   apo_motion, complex_motion, output_dir)
    
    # 7. 生成NMA动画帧
    print("\n" + "="*60)
    print("生成NMA动画帧")
    print("="*60)
    apo_frames_dir = generate_nma_frames(apo_anm, apo_protein, "apo", n_frames=20)
    complex_frames_dir = generate_nma_frames(complex_anm, complex_protein, "complex", n_frames=20)
    
    # 8. 生成汇总表格
    summary_file = output_dir / "nma_summary_table.txt"
    generate_summary_table(apo_eigenvalues, complex_eigenvalues, apo_gate, complex_gate, summary_file)
    
    # 9. 输出eigenvalues到文件
    with open(output_dir / "eigenvalues.txt", 'w') as f:
        f.write("Apo前6个Eigenvalues:\n")
        for i, val in enumerate(apo_eigenvalues[:6]):
            f.write(f"Mode {i+1}: {val:.6f}\n")
        f.write("\nComplex前6个Eigenvalues:\n")
        for i, val in enumerate(complex_eigenvalues[:6]):
            f.write(f"Mode {i+1}: {val:.6f}\n")
    
    print("\n" + "="*80)
    print("✅ 分析完成！")
    print("="*80)
    print(f"\n结果保存在: {output_dir}/")
    print(f"  - eigenvalues.txt: 前6个eigenvalues")
    print(f"  - eigenvalues_comparison.png: Eigenvalues对比图")
    print(f"  - rms_displacement_comparison.png: RMS Displacement对比图")
    print(f"  - nma_summary_table.txt: 数值化汇总表格")
    print(f"  - nma_frames_apo/: Apo的NMA动画帧")
    print(f"  - nma_frames_complex/: Complex的NMA动画帧")
    print("\n⚠️  注意: 请根据实际PDB结构调整关键残基和gate区域的定义！")

if __name__ == "__main__":
    main()

