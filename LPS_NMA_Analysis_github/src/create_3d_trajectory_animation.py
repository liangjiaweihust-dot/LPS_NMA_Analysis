#!/usr/bin/env python3
"""
生成动态3D结构轨迹图
展示不同正常模式下的结构变化轨迹
"""

import prody
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from pathlib import Path
import sys

# Nature风格设置
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'figure.dpi': 100,
    'savefig.dpi': 100,
})

# Nature风格配色
NATURE_COLORS = {
    'apo': '#2E86AB',      # 深蓝色
    'complex': '#A23B72',  # 深紫红色
    'gate1': '#F18F01',    # 橙色
    'gate2': '#C73E1D',    # 红色
    'protein': '#95A5A6',  # 灰色
    'trajectory': '#3498DB', # 蓝色轨迹
    'text': '#2C3E50',
}

def load_structures_and_nma():
    """加载结构并执行NMA分析"""
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
    apo_anm.calcModes(n_modes=6)
    
    print("  执行Complex NMA...")
    complex_anm = prody.ANM('Complex_ANM')
    complex_anm.buildHessian(complex_protein)
    complex_anm.calcModes(n_modes=6)
    
    return apo_protein, complex_protein, apo_anm, complex_anm

def get_backbone_coords(protein, sample=5):
    """获取蛋白质骨架坐标（降采样）"""
    backbone = protein.select('backbone')
    if backbone:
        coords = backbone.getCoords()
        return coords[::sample]
    return protein.getCoords()[::sample*2]

def get_gate_regions(protein):
    """获取gate区域"""
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
        
        gate1_residues = beta1_residues[:10]
        gate2_residues = beta26_residues[-10:]
        
        gate1_atoms = None
        for resnum in gate1_residues:
            res_atoms = protein.select(f'resnum {resnum}')
            if res_atoms is not None:
                if gate1_atoms is None:
                    gate1_atoms = res_atoms
                else:
                    gate1_atoms = gate1_atoms + res_atoms
        
        gate2_atoms = None
        for resnum in gate2_residues:
            res_atoms = protein.select(f'resnum {resnum}')
            if res_atoms is not None:
                if gate2_atoms is None:
                    gate2_atoms = res_atoms
                else:
                    gate2_atoms = gate2_atoms + res_atoms
        
        return gate1_atoms, gate2_atoms
    
    return None, None

def apply_mode_displacement(protein, anm, mode_idx, amplitude):
    """应用正常模式的位移到结构"""
    if mode_idx >= len(anm):
        return None
    
    mode = anm[mode_idx]
    displacements = mode.getArrayNx3()
    
    # 获取原始坐标
    original_coords = protein.getCoords()
    
    # 应用位移
    new_coords = original_coords + displacements * amplitude
    
    return new_coords

def create_3d_trajectory_animation(apo_protein, complex_protein, apo_anm, complex_anm, output_file):
    """创建3D结构轨迹动画"""
    print("\n创建3D结构轨迹动画...")
    
    # 获取骨架坐标
    apo_backbone = get_backbone_coords(apo_protein, sample=8)
    complex_backbone = get_backbone_coords(complex_protein, sample=8)
    
    # 获取gate区域
    apo_gate1, apo_gate2 = get_gate_regions(apo_protein)
    complex_gate1, complex_gate2 = get_gate_regions(complex_protein)
    
    # 创建图形
    fig = plt.figure(figsize=(20, 10))
    
    # 创建2个子图：Apo和Complex
    ax1 = fig.add_subplot(121, projection='3d')  # Apo
    ax2 = fig.add_subplot(122, projection='3d')  # Complex
    
    # 存储轨迹点
    apo_trajectory_points = []
    complex_trajectory_points = []
    apo_gate1_trajectory = []
    apo_gate2_trajectory = []
    complex_gate1_trajectory = []
    complex_gate2_trajectory = []
    
    # 设置坐标范围
    def set_axes_equal(ax, coords):
        if len(coords) == 0:
            return
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
        
        # 计算当前帧的参数
        # 循环显示不同模式：每40帧一个模式
        mode_idx = (frame // 40) % 6
        # 在当前模式内，振幅从-5到5变化
        amplitude = 5 * np.sin((frame % 40) * np.pi / 20)
        
        # ========== Apo结构轨迹 ==========
        # 应用模式位移
        apo_new_coords = apply_mode_displacement(apo_protein, apo_anm, mode_idx, amplitude)
        if apo_new_coords is not None:
            # 获取当前帧的骨架坐标
            backbone_indices = np.arange(0, len(apo_new_coords), 8)
            current_backbone = apo_new_coords[backbone_indices]
            
            # 绘制当前结构（半透明）
            ax1.scatter(current_backbone[:, 0], current_backbone[:, 1], current_backbone[:, 2],
                       c=NATURE_COLORS['apo'], alpha=0.3, s=3, label='Current structure')
            
            # 计算gate中心
            if apo_gate1 is not None and apo_gate2 is not None:
                gate1_indices = apo_gate1.getIndices()
                gate2_indices = apo_gate2.getIndices()
                
                # 找到gate原子在protein中的位置
                protein_indices = apo_protein.getIndices()
                gate1_pos_in_protein = []
                gate2_pos_in_protein = []
                
                for idx in gate1_indices:
                    if idx in protein_indices:
                        pos = np.where(protein_indices == idx)[0]
                        if len(pos) > 0:
                            gate1_pos_in_protein.append(pos[0])
                
                for idx in gate2_indices:
                    if idx in protein_indices:
                        pos = np.where(protein_indices == idx)[0]
                        if len(pos) > 0:
                            gate2_pos_in_protein.append(pos[0])
                
                if gate1_pos_in_protein and gate2_pos_in_protein:
                    gate1_coords = apo_new_coords[gate1_pos_in_protein]
                    gate2_coords = apo_new_coords[gate2_pos_in_protein]
                    
                    gate1_center = gate1_coords.mean(axis=0)
                    gate2_center = gate2_coords.mean(axis=0)
                    
                    # 保存轨迹点
                    apo_gate1_trajectory.append(gate1_center.copy())
                    apo_gate2_trajectory.append(gate2_center.copy())
                    
                    # 只保留最近50个点（轨迹）
                    if len(apo_gate1_trajectory) > 50:
                        apo_gate1_trajectory.pop(0)
                        apo_gate2_trajectory.pop(0)
                    
                    # 绘制gate中心
                    ax1.scatter([gate1_center[0]], [gate1_center[1]], [gate1_center[2]],
                               c=NATURE_COLORS['gate1'], s=400, marker='o',
                               edgecolors='black', linewidths=2, label='Gate 1')
                    ax1.scatter([gate2_center[0]], [gate2_center[1]], [gate2_center[2]],
                               c=NATURE_COLORS['gate2'], s=400, marker='s',
                               edgecolors='black', linewidths=2, label='Gate 2')
                    
                    # 绘制gate距离线
                    ax1.plot([gate1_center[0], gate2_center[0]],
                            [gate1_center[1], gate2_center[1]],
                            [gate1_center[2], gate2_center[2]],
                            'k-', linewidth=2, alpha=0.6)
                    
                    # 绘制轨迹
                    if len(apo_gate1_trajectory) > 1:
                        traj1 = np.array(apo_gate1_trajectory)
                        traj2 = np.array(apo_gate2_trajectory)
                        
                        # Gate 1轨迹
                        ax1.plot(traj1[:, 0], traj1[:, 1], traj1[:, 2],
                                color=NATURE_COLORS['gate1'], linewidth=2, alpha=0.6,
                                linestyle='--', label='Gate 1 trajectory')
                        
                        # Gate 2轨迹
                        ax1.plot(traj2[:, 0], traj2[:, 1], traj2[:, 2],
                                color=NATURE_COLORS['gate2'], linewidth=2, alpha=0.6,
                                linestyle='--', label='Gate 2 trajectory')
                    
                    # 计算距离
                    distance = np.linalg.norm(gate2_center - gate1_center)
                    
                    # 显示信息
                    info_text = f'Apo Structure\n'
                    info_text += f'Mode {mode_idx+1}\n'
                    info_text += f'Amplitude: {amplitude:.2f}\n'
                    info_text += f'Gate Distance: {distance:.2f} Å'
                    
                    ax1.text2D(0.05, 0.95, info_text,
                             transform=ax1.transAxes, fontsize=11, fontweight='bold',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        
        ax1.set_title('Apo Structure Trajectory', fontweight='bold', fontsize=13,
                     color=NATURE_COLORS['text'], pad=15)
        set_axes_equal(ax1, apo_backbone)
        ax1.set_xlabel('X (Å)', fontsize=10)
        ax1.set_ylabel('Y (Å)', fontsize=10)
        ax1.set_zlabel('Z (Å)', fontsize=10)
        if frame == 0:  # 只在第一帧显示图例
            ax1.legend(loc='upper left', fontsize=9, framealpha=0.9)
        
        # ========== Complex结构轨迹 ==========
        complex_new_coords = apply_mode_displacement(complex_protein, complex_anm, mode_idx, amplitude)
        if complex_new_coords is not None:
            backbone_indices = np.arange(0, len(complex_new_coords), 8)
            current_backbone = complex_new_coords[backbone_indices]
            
            ax2.scatter(current_backbone[:, 0], current_backbone[:, 1], current_backbone[:, 2],
                       c=NATURE_COLORS['complex'], alpha=0.3, s=3)
            
            if complex_gate1 is not None and complex_gate2 is not None:
                gate1_indices = complex_gate1.getIndices()
                gate2_indices = complex_gate2.getIndices()
                
                protein_indices = complex_protein.getIndices()
                gate1_pos_in_protein = []
                gate2_pos_in_protein = []
                
                for idx in gate1_indices:
                    if idx in protein_indices:
                        pos = np.where(protein_indices == idx)[0]
                        if len(pos) > 0:
                            gate1_pos_in_protein.append(pos[0])
                
                for idx in gate2_indices:
                    if idx in protein_indices:
                        pos = np.where(protein_indices == idx)[0]
                        if len(pos) > 0:
                            gate2_pos_in_protein.append(pos[0])
                
                if gate1_pos_in_protein and gate2_pos_in_protein:
                    gate1_coords = complex_new_coords[gate1_pos_in_protein]
                    gate2_coords = complex_new_coords[gate2_pos_in_protein]
                    
                    gate1_center = gate1_coords.mean(axis=0)
                    gate2_center = gate2_coords.mean(axis=0)
                    
                    complex_gate1_trajectory.append(gate1_center.copy())
                    complex_gate2_trajectory.append(gate2_center.copy())
                    
                    if len(complex_gate1_trajectory) > 50:
                        complex_gate1_trajectory.pop(0)
                        complex_gate2_trajectory.pop(0)
                    
                    ax2.scatter([gate1_center[0]], [gate1_center[1]], [gate1_center[2]],
                               c=NATURE_COLORS['gate1'], s=400, marker='o',
                               edgecolors='black', linewidths=2)
                    ax2.scatter([gate2_center[0]], [gate2_center[1]], [gate2_center[2]],
                               c=NATURE_COLORS['gate2'], s=400, marker='s',
                               edgecolors='black', linewidths=2)
                    
                    ax2.plot([gate1_center[0], gate2_center[0]],
                            [gate1_center[1], gate2_center[1]],
                            [gate1_center[2], gate2_center[2]],
                            'k-', linewidth=2, alpha=0.6)
                    
                    if len(complex_gate1_trajectory) > 1:
                        traj1 = np.array(complex_gate1_trajectory)
                        traj2 = np.array(complex_gate2_trajectory)
                        
                        ax2.plot(traj1[:, 0], traj1[:, 1], traj1[:, 2],
                                color=NATURE_COLORS['gate1'], linewidth=2, alpha=0.6,
                                linestyle='--')
                        ax2.plot(traj2[:, 0], traj2[:, 1], traj2[:, 2],
                                color=NATURE_COLORS['gate2'], linewidth=2, alpha=0.6,
                                linestyle='--')
                    
                    distance = np.linalg.norm(gate2_center - gate1_center)
                    
                    info_text = f'Complex Structure\n'
                    info_text += f'Mode {mode_idx+1}\n'
                    info_text += f'Amplitude: {amplitude:.2f}\n'
                    info_text += f'Gate Distance: {distance:.2f} Å'
                    
                    ax2.text2D(0.05, 0.95, info_text,
                             transform=ax2.transAxes, fontsize=11, fontweight='bold',
                             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        
        ax2.set_title('Complex Structure Trajectory', fontweight='bold', fontsize=13,
                     color=NATURE_COLORS['text'], pad=15)
        set_axes_equal(ax2, complex_backbone)
        ax2.set_xlabel('X (Å)', fontsize=10)
        ax2.set_ylabel('Y (Å)', fontsize=10)
        ax2.set_zlabel('Z (Å)', fontsize=10)
        
        # 添加总标题
        fig.suptitle(f'3D Structure Trajectory Animation - Frame {frame+1}/240', 
                    fontsize=15, fontweight='bold', y=0.98)
        
        return ax1, ax2
    
    # 创建动画
    print("  生成动画帧（这可能需要一些时间）...")
    anim = animation.FuncAnimation(fig, animate, frames=240, interval=80, blit=False)
    
    # 保存动画
    print("  保存3D轨迹动画...")
    anim.save(output_file, writer='pillow', fps=12, dpi=80)
    print(f"  ✅ 3D轨迹动画已保存: {output_file}")
    
    plt.close()

def create_mode_comparison_trajectory(apo_protein, complex_protein, apo_anm, complex_anm, output_file):
    """创建不同模式的轨迹对比图"""
    print("\n创建不同模式的轨迹对比图...")
    
    # 获取gate区域
    apo_gate1, apo_gate2 = get_gate_regions(apo_protein)
    complex_gate1, complex_gate2 = get_gate_regions(complex_protein)
    
    fig = plt.figure(figsize=(20, 12))
    
    # 创建6个子图，每个模式一个
    axes = []
    for i in range(6):
        ax = fig.add_subplot(2, 3, i+1, projection='3d')
        axes.append(ax)
    
    # 为每个模式生成轨迹
    n_points = 50
    amplitudes = np.linspace(-5, 5, n_points)
    
    for mode_idx in range(6):
        ax = axes[mode_idx]
        
        # Apo轨迹
        apo_gate1_traj = []
        apo_gate2_traj = []
        complex_gate1_traj = []
        complex_gate2_traj = []
        
        for amp in amplitudes:
            # Apo
            apo_new_coords = apply_mode_displacement(apo_protein, apo_anm, mode_idx, amp)
            if apo_new_coords is not None and apo_gate1 is not None and apo_gate2 is not None:
                gate1_indices = apo_gate1.getIndices()
                gate2_indices = apo_gate2.getIndices()
                protein_indices = apo_protein.getIndices()
                
                gate1_pos = [np.where(protein_indices == idx)[0][0] 
                            for idx in gate1_indices if idx in protein_indices]
                gate2_pos = [np.where(protein_indices == idx)[0][0] 
                            for idx in gate2_indices if idx in protein_indices]
                
                if gate1_pos and gate2_pos:
                    gate1_center = apo_new_coords[gate1_pos].mean(axis=0)
                    gate2_center = apo_new_coords[gate2_pos].mean(axis=0)
                    apo_gate1_traj.append(gate1_center)
                    apo_gate2_traj.append(gate2_center)
            
            # Complex
            complex_new_coords = apply_mode_displacement(complex_protein, complex_anm, mode_idx, amp)
            if complex_new_coords is not None and complex_gate1 is not None and complex_gate2 is not None:
                gate1_indices = complex_gate1.getIndices()
                gate2_indices = complex_gate2.getIndices()
                protein_indices = complex_protein.getIndices()
                
                gate1_pos = [np.where(protein_indices == idx)[0][0] 
                            for idx in gate1_indices if idx in protein_indices]
                gate2_pos = [np.where(protein_indices == idx)[0][0] 
                            for idx in gate2_indices if idx in protein_indices]
                
                if gate1_pos and gate2_pos:
                    gate1_center = complex_new_coords[gate1_pos].mean(axis=0)
                    gate2_center = complex_new_coords[gate2_pos].mean(axis=0)
                    complex_gate1_traj.append(gate1_center)
                    complex_gate2_traj.append(gate2_center)
        
        # 绘制轨迹
        if apo_gate1_traj and apo_gate2_traj:
            apo_traj1 = np.array(apo_gate1_traj)
            apo_traj2 = np.array(apo_gate2_traj)
            
            ax.plot(apo_traj1[:, 0], apo_traj1[:, 1], apo_traj1[:, 2],
                   color=NATURE_COLORS['apo'], linewidth=2.5, alpha=0.8,
                   label='Apo Gate 1', linestyle='-')
            ax.plot(apo_traj2[:, 0], apo_traj2[:, 1], apo_traj2[:, 2],
                   color=NATURE_COLORS['apo'], linewidth=2.5, alpha=0.8,
                   label='Apo Gate 2', linestyle='--')
            
            # 起点和终点
            ax.scatter([apo_traj1[0, 0]], [apo_traj1[0, 1]], [apo_traj1[0, 2]],
                      c=NATURE_COLORS['apo'], s=200, marker='o', edgecolors='black', linewidths=2)
            ax.scatter([apo_traj1[-1, 0]], [apo_traj1[-1, 1]], [apo_traj1[-1, 2]],
                      c=NATURE_COLORS['apo'], s=200, marker='s', edgecolors='black', linewidths=2)
        
        if complex_gate1_traj and complex_gate2_traj:
            complex_traj1 = np.array(complex_gate1_traj)
            complex_traj2 = np.array(complex_gate2_traj)
            
            ax.plot(complex_traj1[:, 0], complex_traj1[:, 1], complex_traj1[:, 2],
                   color=NATURE_COLORS['complex'], linewidth=2.5, alpha=0.8,
                   label='Complex Gate 1', linestyle='-')
            ax.plot(complex_traj2[:, 0], complex_traj2[:, 1], complex_traj2[:, 2],
                   color=NATURE_COLORS['complex'], linewidth=2.5, alpha=0.8,
                   label='Complex Gate 2', linestyle='--')
            
            ax.scatter([complex_traj1[0, 0]], [complex_traj1[0, 1]], [complex_traj1[0, 2]],
                      c=NATURE_COLORS['complex'], s=200, marker='o', edgecolors='black', linewidths=2)
            ax.scatter([complex_traj1[-1, 0]], [complex_traj1[-1, 1]], [complex_traj1[-1, 2]],
                      c=NATURE_COLORS['complex'], s=200, marker='s', edgecolors='black', linewidths=2)
        
        ax.set_title(f'Mode {mode_idx+1} Trajectory', fontweight='bold', fontsize=12)
        ax.set_xlabel('X (Å)', fontsize=9)
        ax.set_ylabel('Y (Å)', fontsize=9)
        ax.set_zlabel('Z (Å)', fontsize=9)
        if mode_idx == 0:
            ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
    
    plt.suptitle('Gate Region Trajectories for Different Normal Modes', 
                fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, facecolor='white')
    plt.close()
    print(f"  ✅ 模式对比轨迹图已保存: {output_file}")

def main():
    """主函数"""
    print("="*80)
    print("3D结构轨迹动画生成")
    print("="*80)
    
    # 加载数据
    apo_protein, complex_protein, apo_anm, complex_anm = load_structures_and_nma()
    
    # 创建输出目录
    output_dir = Path("nma_results")
    output_dir.mkdir(exist_ok=True)
    
    # 创建3D轨迹动画
    create_3d_trajectory_animation(apo_protein, complex_protein, apo_anm, complex_anm,
                                  output_dir / "structure_trajectory_3d.gif")
    
    # 创建模式对比轨迹图
    create_mode_comparison_trajectory(apo_protein, complex_protein, apo_anm, complex_anm,
                                     output_dir / "mode_trajectory_comparison.png")
    
    print("\n" + "="*80)
    print("✅ 所有轨迹可视化已完成！")
    print("="*80)
    print(f"\n生成的文件:")
    print(f"  - {output_dir}/structure_trajectory_3d.gif (动态3D结构轨迹动画)")
    print(f"  - {output_dir}/mode_trajectory_comparison.png (不同模式的轨迹对比图)")
    print("\n图例说明:")
    print("  🟠 橙色球体 (○) = Gate 1 (β1区域)")
    print("  🔴 红色方块 (□) = Gate 2 (β26区域)")
    print("  ➖ 虚线 = Gate轨迹路径")
    print("  📍 起点 = 初始位置")
    print("  📍 终点 = 最大位移位置")

if __name__ == "__main__":
    main()



