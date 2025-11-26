#!/usr/bin/env python3
"""
高级Gate Openness可视化脚本
生成Nature风格的图表和3D动画
"""

import prody
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from pathlib import Path
import sys

# Nature风格设置
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 1.2,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.minor.width': 0.8,
    'ytick.minor.width': 0.8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'legend.frameon': False
})

# Nature风格配色方案
NATURE_COLORS = {
    'apo': '#2E86AB',      # 深蓝色
    'complex': '#A23B72',  # 深紫红色
    'background': '#FFFFFF',
    'grid': '#E8E8E8',
    'text': '#2C3E50',
    'accent1': '#F18F01',  # 橙色
    'accent2': '#C73E1D',  # 红色
}

def load_nma_results():
    """从之前的分析结果中加载数据"""
    # 重新运行NMA分析获取gate openness数据
    print("加载结构并执行NMA分析...")
    
    apo_pdb = "no-peptide.pdb"
    complex_pdb = "complex.pdb"
    
    # 加载结构
    apo_structure = prody.parsePDB(apo_pdb)
    complex_structure = prody.parsePDB(complex_pdb)
    
    apo_protein = apo_structure.select('protein')
    complex_protein = complex_structure.select('protein')
    
    # 执行NMA
    print("  执行Apo NMA...")
    apo_anm = prody.ANM('Apo_ANM')
    apo_anm.buildHessian(apo_protein)
    apo_anm.calcModes(n_modes=20)
    
    print("  执行Complex NMA...")
    complex_anm = prody.ANM('Complex_ANM')
    complex_anm.buildHessian(complex_protein)
    complex_anm.calcModes(n_modes=20)
    
    # 识别gate区域（使用与主脚本相同的逻辑）
    def get_gate_residues(protein):
        chains = sorted(set(protein.getChids()))
        chain_info = {}
        for chid in chains:
            chain = protein.select(f'chain {chid}')
            if chain:
                resnums = sorted(set(chain.getResnums()))
                chain_info[chid] = resnums
        
        if chain_info:
            longest_chain = max(chain_info.items(), key=lambda x: len(x[1]))
            main_chain = longest_chain[0]
            main_resnums = chain_info[main_chain]
            
            beta1_residues = main_resnums[:len(main_resnums)//4][:20]
            beta26_residues = main_resnums[3*len(main_resnums)//4:][-20:]
            
            gate1 = beta1_residues[:10]
            gate2 = beta26_residues[-10:]
            return gate1, gate2
        return [], []
    
    apo_gate1, apo_gate2 = get_gate_residues(apo_protein)
    complex_gate1, complex_gate2 = get_gate_residues(complex_protein)
    
    # 计算gate openness
    def calculate_gate_openness_detailed(anm, protein, gate_residues_1, gate_residues_2):
        """计算详细的gate openness数据"""
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
            return None
        
        gate1_center = gate1_atoms.getCoords().mean(axis=0)
        gate2_center = gate2_atoms.getCoords().mean(axis=0)
        initial_distance = np.linalg.norm(gate2_center - gate1_center)
        
        protein_indices = protein.getIndices()
        gate1_selection_indices = gate1_atoms.getIndices()
        gate2_selection_indices = gate2_atoms.getIndices()
        
        gate_openness_per_mode = []
        gate1_coords_per_mode = []
        gate2_coords_per_mode = []
        
        for mode_idx in range(min(6, len(anm))):
            mode = anm[mode_idx]
            displacements = mode.getArrayNx3()
            
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
            
            gate1_valid = gate1_indices_in_protein[gate1_indices_in_protein < len(displacements)]
            gate2_valid = gate2_indices_in_protein[gate2_indices_in_protein < len(displacements)]
            
            if len(gate1_valid) == 0 or len(gate2_valid) == 0:
                continue
            
            gate1_disp = displacements[gate1_valid].mean(axis=0)
            gate2_disp = displacements[gate2_valid].mean(axis=0)
            
            gate_vector = gate2_center - gate1_center
            if np.linalg.norm(gate_vector) < 1e-6:
                continue
            gate_vector_norm = gate_vector / np.linalg.norm(gate_vector)
            
            gate1_proj = np.dot(gate1_disp, -gate_vector_norm)
            gate2_proj = np.dot(gate2_disp, gate_vector_norm)
            openness = gate1_proj + gate2_proj
            
            gate_openness_per_mode.append(openness)
            gate1_coords_per_mode.append(gate1_center + gate1_disp * 5)  # 放大5倍用于可视化
            gate2_coords_per_mode.append(gate2_center + gate2_disp * 5)
        
        return {
            'initial_distance': initial_distance,
            'gate1_center': gate1_center,
            'gate2_center': gate2_center,
            'openness_per_mode': gate_openness_per_mode,
            'gate1_coords_per_mode': gate1_coords_per_mode,
            'gate2_coords_per_mode': gate2_coords_per_mode,
            'average_openness': np.mean(gate_openness_per_mode) if gate_openness_per_mode else 0,
            'max_openness': np.max(gate_openness_per_mode) if gate_openness_per_mode else 0
        }
    
    print("  计算Gate Openness...")
    apo_gate_data = calculate_gate_openness_detailed(apo_anm, apo_protein, apo_gate1, apo_gate2)
    complex_gate_data = calculate_gate_openness_detailed(complex_anm, complex_protein, complex_gate1, complex_gate2)
    
    return apo_gate_data, complex_gate_data, apo_protein, complex_protein

def create_advanced_gate_plot(apo_data, complex_data, output_file):
    """创建高级Gate Openness可视化图"""
    print("\n创建高级Gate Openness可视化图...")
    
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3, 
                          left=0.08, right=0.95, top=0.95, bottom=0.08)
    
    # 子图1: Gate Openness对比（柱状图）
    ax1 = fig.add_subplot(gs[0, 0])
    modes = np.arange(1, len(apo_data['openness_per_mode']) + 1)
    width = 0.35
    x = np.arange(len(modes))
    
    bars1 = ax1.bar(x - width/2, apo_data['openness_per_mode'], width, 
                    label='Apo', color=NATURE_COLORS['apo'], alpha=0.85, 
                    edgecolor='white', linewidth=1.5)
    bars2 = ax1.bar(x + width/2, complex_data['openness_per_mode'], width, 
                    label='Complex', color=NATURE_COLORS['complex'], alpha=0.85,
                    edgecolor='white', linewidth=1.5)
    
    # 添加数值标签
    for i, (bar1, bar2) in enumerate(zip(bars1, bars2)):
        height1 = bar1.get_height()
        height2 = bar2.get_height()
        if abs(height1) > 0.001:
            ax1.text(bar1.get_x() + bar1.get_width()/2., height1,
                    f'{height1:.3f}', ha='center', va='bottom' if height1 > 0 else 'top',
                    fontsize=8, color=NATURE_COLORS['apo'])
        if abs(height2) > 0.001:
            ax1.text(bar2.get_x() + bar2.get_width()/2., height2,
                    f'{height2:.3f}', ha='center', va='bottom' if height2 > 0 else 'top',
                    fontsize=8, color=NATURE_COLORS['complex'])
    
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.8, alpha=0.3)
    ax1.set_xlabel('Mode Number', fontweight='bold', color=NATURE_COLORS['text'])
    ax1.set_ylabel('Gate Openness (Å)', fontweight='bold', color=NATURE_COLORS['text'])
    ax1.set_title('Gate Openness per Mode', fontweight='bold', fontsize=13, 
                  color=NATURE_COLORS['text'], pad=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels(modes)
    ax1.legend(loc='upper right', frameon=True, fancybox=True, shadow=True, 
              framealpha=0.9, edgecolor='none')
    ax1.grid(True, alpha=0.2, linestyle='--', linewidth=0.8)
    
    # 子图2: 平均Gate Openness对比
    ax2 = fig.add_subplot(gs[0, 1])
    categories = ['Apo', 'Complex']
    values = [apo_data['average_openness'], complex_data['average_openness']]
    colors = [NATURE_COLORS['apo'], NATURE_COLORS['complex']]
    
    bars = ax2.bar(categories, values, color=colors, alpha=0.85, 
                   edgecolor='white', linewidth=2, width=0.6)
    
    # 添加数值标签和趋势箭头
    for i, (bar, val) in enumerate(zip(bars, values)):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.4f}', ha='center', 
                va='bottom' if height > 0 else 'top',
                fontsize=11, fontweight='bold', color=colors[i])
        
        # 添加趋势箭头
        if i == 0 and height > 0:
            ax2.annotate('', xy=(bar.get_x() + bar.get_width()/2, height + 0.001),
                        xytext=(bar.get_x() + bar.get_width()/2, height),
                        arrowprops=dict(arrowstyle='->', lw=2, color=NATURE_COLORS['accent1']))
            ax2.text(bar.get_x() + bar.get_width()/2, height + 0.002,
                    'Open', ha='center', fontsize=9, color=NATURE_COLORS['accent1'],
                    fontweight='bold')
        elif i == 1 and height < 0:
            ax2.annotate('', xy=(bar.get_x() + bar.get_width()/2, height - 0.001),
                        xytext=(bar.get_x() + bar.get_width()/2, height),
                        arrowprops=dict(arrowstyle='->', lw=2, color=NATURE_COLORS['accent2']))
            ax2.text(bar.get_x() + bar.get_width()/2, height - 0.002,
                    'Closed', ha='center', fontsize=9, color=NATURE_COLORS['accent2'],
                    fontweight='bold')
    
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    ax2.set_ylabel('Average Gate Openness (Å)', fontweight='bold', 
                   color=NATURE_COLORS['text'])
    ax2.set_title('Average Gate Openness Comparison', fontweight='bold', 
                 fontsize=13, color=NATURE_COLORS['text'], pad=10)
    ax2.grid(True, alpha=0.2, linestyle='--', linewidth=0.8, axis='y')
    
    # 子图3: Gate距离变化
    ax3 = fig.add_subplot(gs[1, 0])
    initial_distances = [apo_data['initial_distance'], complex_data['initial_distance']]
    bars = ax3.bar(categories, initial_distances, color=colors, alpha=0.85,
                   edgecolor='white', linewidth=2, width=0.6)
    
    for bar, val in zip(bars, initial_distances):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height,
                f'{val:.2f} Å', ha='center', va='bottom',
                fontsize=11, fontweight='bold', color='white')
    
    ax3.set_ylabel('Initial Gate Distance (Å)', fontweight='bold',
                   color=NATURE_COLORS['text'])
    ax3.set_title('Initial Gate Distance', fontweight='bold',
                 fontsize=13, color=NATURE_COLORS['text'], pad=10)
    ax3.grid(True, alpha=0.2, linestyle='--', linewidth=0.8, axis='y')
    
    # 子图4: 综合趋势图
    ax4 = fig.add_subplot(gs[1, 1])
    
    # 创建综合指标（归一化）
    apo_norm = (apo_data['average_openness'] - min(apo_data['average_openness'], 
              complex_data['average_openness'])) / (max(apo_data['average_openness'],
              complex_data['average_openness']) - min(apo_data['average_openness'],
              complex_data['average_openness']) + 1e-6) if abs(apo_data['average_openness'] - 
              complex_data['average_openness']) > 1e-6 else 0.5
    
    complex_norm = 1 - apo_norm
    
    # 雷达图风格
    angles = np.linspace(0, 2*np.pi, 4, endpoint=False).tolist()
    angles += angles[:1]  # 闭合
    
    apo_values = [apo_norm, 1 - abs(apo_data['average_openness']), 
                  apo_data['initial_distance']/100, 0.5]
    complex_values = [complex_norm, 1 - abs(complex_data['average_openness']),
                     complex_data['initial_distance']/100, 0.5]
    
    apo_values += apo_values[:1]
    complex_values += complex_values[:1]
    
    ax4.plot(angles, apo_values, 'o-', linewidth=2.5, label='Apo',
            color=NATURE_COLORS['apo'], markersize=8, alpha=0.8)
    ax4.fill(angles, apo_values, alpha=0.25, color=NATURE_COLORS['apo'])
    
    ax4.plot(angles, complex_values, 's-', linewidth=2.5, label='Complex',
            color=NATURE_COLORS['complex'], markersize=8, alpha=0.8)
    ax4.fill(angles, complex_values, alpha=0.25, color=NATURE_COLORS['complex'])
    
    ax4.set_xticks(angles[:-1])
    ax4.set_xticklabels(['Openness', 'Stability', 'Distance', 'Rigidity'],
                       fontsize=9)
    ax4.set_ylim(0, 1)
    ax4.set_title('Gate Properties Comparison', fontweight='bold',
                 fontsize=13, color=NATURE_COLORS['text'], pad=15)
    ax4.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
    ax4.grid(True, alpha=0.3, linestyle='--')
    
    # 添加总标题
    fig.suptitle('Gate Openness Analysis: Apo vs Complex', 
                fontsize=16, fontweight='bold', 
                color=NATURE_COLORS['text'], y=0.98)
    
    plt.savefig(output_file, dpi=300, facecolor='white', edgecolor='none')
    plt.close()
    print(f"  ✅ 高级Gate Openness图已保存: {output_file}")

def create_3d_gate_animation(apo_data, complex_data, apo_protein, complex_protein, output_file):
    """创建3D Gate Openness动画"""
    print("\n创建3D Gate Openness动画...")
    
    fig = plt.figure(figsize=(16, 8))
    
    # 创建两个3D子图
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')
    
    # 获取蛋白质骨架坐标（用于背景）
    def get_backbone_coords(protein):
        backbone = protein.select('backbone')
        if backbone:
            return backbone.getCoords()
        return protein.getCoords()[::10]  # 降采样
    
    apo_backbone = get_backbone_coords(apo_protein)
    complex_backbone = get_backbone_coords(complex_protein)
    
    # 设置坐标范围
    def set_axes_equal(ax, coords):
        max_range = np.array([coords[:, i].max() - coords[:, i].min() 
                             for i in range(3)]).max() / 2.0
        mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
        mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
        mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    # 动画函数
    def animate(frame):
        ax1.clear()
        ax2.clear()
        
        # 计算当前帧的振幅（-5到5，循环）
        amplitude = 5 * np.sin(frame * np.pi / 30)
        
        # Apo结构
        ax1.scatter(apo_backbone[:, 0], apo_backbone[:, 1], apo_backbone[:, 2],
                   c=NATURE_COLORS['apo'], alpha=0.1, s=1)
        
        # Gate区域
        gate1_center = apo_data['gate1_center']
        gate2_center = apo_data['gate2_center']
        
        # 计算当前模式的位移
        mode_idx = frame % len(apo_data['openness_per_mode'])
        if mode_idx < len(apo_data['gate1_coords_per_mode']):
            gate1_pos = gate1_center + (apo_data['gate1_coords_per_mode'][mode_idx] - 
                                       gate1_center) * amplitude / 5
            gate2_pos = gate2_center + (apo_data['gate2_coords_per_mode'][mode_idx] - 
                                       gate2_center) * amplitude / 5
        else:
            gate1_pos = gate1_center
            gate2_pos = gate2_center
        
        # 绘制gate区域
        ax1.scatter([gate1_pos[0]], [gate1_pos[1]], [gate1_pos[2]],
                   c=NATURE_COLORS['accent1'], s=200, marker='o', 
                   edgecolors='white', linewidths=2, label='Gate 1')
        ax1.scatter([gate2_pos[0]], [gate2_pos[1]], [gate2_pos[2]],
                   c=NATURE_COLORS['accent2'], s=200, marker='s',
                   edgecolors='white', linewidths=2, label='Gate 2')
        
        # 绘制连接线
        ax1.plot([gate1_pos[0], gate2_pos[0]], 
                [gate1_pos[1], gate2_pos[1]],
                [gate1_pos[2], gate2_pos[2]],
                'k--', linewidth=2, alpha=0.6)
        
        # 计算距离
        distance = np.linalg.norm(gate2_pos - gate1_pos)
        openness = apo_data['openness_per_mode'][mode_idx] * amplitude / 5
        
        ax1.set_title(f'Apo - Mode {mode_idx+1}\nDistance: {distance:.2f} Å, '
                     f'Openness: {openness:.3f} Å',
                     fontweight='bold', fontsize=11, color=NATURE_COLORS['text'])
        set_axes_equal(ax1, apo_backbone)
        ax1.set_xlabel('X (Å)', fontsize=9)
        ax1.set_ylabel('Y (Å)', fontsize=9)
        ax1.set_zlabel('Z (Å)', fontsize=9)
        
        # Complex结构
        ax2.scatter(complex_backbone[:, 0], complex_backbone[:, 1], complex_backbone[:, 2],
                    c=NATURE_COLORS['complex'], alpha=0.1, s=1)
        
        gate1_center_c = complex_data['gate1_center']
        gate2_center_c = complex_data['gate2_center']
        
        if mode_idx < len(complex_data['gate1_coords_per_mode']):
            gate1_pos_c = gate1_center_c + (complex_data['gate1_coords_per_mode'][mode_idx] - 
                                           gate1_center_c) * amplitude / 5
            gate2_pos_c = gate2_center_c + (complex_data['gate2_coords_per_mode'][mode_idx] - 
                                            gate2_center_c) * amplitude / 5
        else:
            gate1_pos_c = gate1_center_c
            gate2_pos_c = gate2_center_c
        
        ax2.scatter([gate1_pos_c[0]], [gate1_pos_c[1]], [gate1_pos_c[2]],
                   c=NATURE_COLORS['accent1'], s=200, marker='o',
                   edgecolors='white', linewidths=2)
        ax2.scatter([gate2_pos_c[0]], [gate2_pos_c[1]], [gate2_pos_c[2]],
                   c=NATURE_COLORS['accent2'], s=200, marker='s',
                   edgecolors='white', linewidths=2)
        
        ax2.plot([gate1_pos_c[0], gate2_pos_c[0]],
                [gate1_pos_c[1], gate2_pos_c[1]],
                [gate1_pos_c[2], gate2_pos_c[2]],
                'k--', linewidth=2, alpha=0.6)
        
        distance_c = np.linalg.norm(gate2_pos_c - gate1_pos_c)
        openness_c = complex_data['openness_per_mode'][mode_idx] * amplitude / 5
        
        ax2.set_title(f'Complex - Mode {mode_idx+1}\nDistance: {distance_c:.2f} Å, '
                     f'Openness: {openness_c:.3f} Å',
                     fontweight='bold', fontsize=11, color=NATURE_COLORS['text'])
        set_axes_equal(ax2, complex_backbone)
        ax2.set_xlabel('X (Å)', fontsize=9)
        ax2.set_ylabel('Y (Å)', fontsize=9)
        ax2.set_zlabel('Z (Å)', fontsize=9)
        
        return ax1, ax2
    
    # 创建动画
    anim = animation.FuncAnimation(fig, animate, frames=120, interval=100, blit=False)
    
    # 保存动画
    print("  保存3D动画（这可能需要一些时间）...")
    anim.save(output_file, writer='pillow', fps=10, dpi=100)
    print(f"  ✅ 3D动画已保存: {output_file}")
    
    plt.close()

def main():
    """主函数"""
    print("="*80)
    print("高级Gate Openness可视化")
    print("="*80)
    
    # 加载数据
    apo_data, complex_data, apo_protein, complex_protein = load_nma_results()
    
    if apo_data is None or complex_data is None:
        print("错误: 无法计算gate openness数据")
        sys.exit(1)
    
    # 创建输出目录
    output_dir = Path("nma_results")
    output_dir.mkdir(exist_ok=True)
    
    # 创建高级图表
    create_advanced_gate_plot(apo_data, complex_data, 
                            output_dir / "gate_openness_advanced.png")
    
    # 创建3D动画
    create_3d_gate_animation(apo_data, complex_data, apo_protein, complex_protein,
                            output_dir / "gate_openness_3d_animation.gif")
    
    print("\n" + "="*80)
    print("✅ 所有可视化已完成！")
    print("="*80)
    print(f"\n生成的文件:")
    print(f"  - {output_dir}/gate_openness_advanced.png")
    print(f"  - {output_dir}/gate_openness_3d_animation.gif")

if __name__ == "__main__":
    main()



