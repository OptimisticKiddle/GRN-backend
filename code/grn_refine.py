import torch
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from sklearn.preprocessing import normalize
import scipy.io.mmio  # 引入读取mm格式文件的模块
import pandas as pd

import math
import numpy as np
from scipy.stats import pearsonr

# 用于快速计算各类指标 mse、pearson corr、R²、F范数等
from torchmetrics.functional.regression import spearman_corrcoef,r2_score,mean_absolute_percentage_error
from torchmetrics.functional.regression import pearson_corrcoef,cosine_similarity
from torchmetrics.functional.pairwise import pairwise_euclidean_distance,pairwise_manhattan_distance

def load_mm_data(filename):
    return scipy.sparse.coo_matrix(scipy.io.mmread(filename))
def grnRefine(gse_id:str,gsm_id:str):
	basePath = f'../static/GSE{gse_id}/GSM{gsm_id}'
    # 检查GPU 可用
	device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")	
	# 各文件路径
	# wt、ko 的gene by cluster格式的GAM
	wt_csv_path = basePath + '/refine/wt_avg_expr_by_cluster.csv'
	ko_csv_path = basePath + '/refine/ko_avg_expr_by_cluster.csv'
	GRN_path = basePath + '/GRN.mm'
	gene_name = ["Tet2"] # 涉及到的敲除gene symbol 
	num_iterations = 2  # 模拟扰动效应 传播次数 默认2
	# 加载数据:
	# 读取两个CSV文件，同时指定header=None表示没有列名
	wt_data = pd.read_csv(wt_csv_path, index_col=0)
	ko_data = pd.read_csv(ko_csv_path, index_col=0)
	# 对WT数据进行转置并重置索引
	wt_data_transposed = wt_data.T
	# 将转置后的WT数据和原始的KO数据转换成PyTorch张量，设置数据类型为float
	wt_mat_t = torch.tensor(wt_data_transposed.values, dtype=torch.float).to(device)
	ko_mat_t = torch.tensor(ko_data.values, dtype=torch.float).t().to(device)
	wt_mat = torch.tensor(wt_data.values, dtype=torch.float).to(device)
	ko_mat = torch.tensor(ko_data.values, dtype=torch.float).to(device)
    # 加载GRN_init矩阵
	GRN_init_coo = load_mm_data(GRN_path)
	# 转换为稠密张量（这里不再添加偏置）
	GRN_init = torch.tensor(GRN_init_coo.todense(), dtype=torch.float).to(device)
	# 定义模型参数，这里直接使用GRN作为可学习参数
	GRN = nn.Parameter(GRN_init.clone().requires_grad_()).to(device)
    #训练：
    # 找到wt_mat_t中 KO 基因列的索引
	gene_index = []
	for i,x in enumerate(gene_name):
		gene_index.append(wt_data.index.get_loc(x)) # 该基因在wt_mat_t的列索引

	# 获取ko gene列的相反数
	init_col = []
	for _,idx in enumerate(gene_index):
		init_col.append(-wt_mat_t[:, idx].to(device))
	# 定义记录损失的列表
	mse_losses = [] # mse
	reg_losses = []
	total_losses = []
	# 定义 平均相关性
	pearson_avg_corrs = []
	spearman_avg_corrs = []
	cos_avg_corrs = []
	# 欧式距离
	euclidean_avg_dists = []
	# 曼哈顿距离
	manhattan_avg_dists = []
	# 误差 
	mae_lists = [] # 平均绝对误差
	RMSE_lists = [] # 均方根误差
	MAPE_lists = [] # 平均绝对百分比误差
	# R square
	r2_score_list = []
	# 定义MAE损失函数
	mae_criterion = nn.L1Loss(reduction='mean')
	# 定义损失函数
	class MatrixRefinementLoss(nn.Module):
		def __init__(self, init_col, wt_mat_t, reg_weight=1e-3):
			super().__init__()
			self.init_col = init_col
			self.wt_mat_t = wt_mat_t
			self.reg_weight = reg_weight

		def forward(self):
			# 模拟扰动 In Silico Gene Perturbation
			# 初始化delta_mat  
			delta_mat = torch.zeros_like(self.wt_mat_t)
			for i,idx in enumerate(gene_index):
				delta_mat[:, idx] = self.init_col[i].clone()

			# 信号传播 Signal Propagation
			# 进行两次循环计算
			for _ in range(num_iterations):
				# 矩阵乘法更新delta_mat
				delta_mat = torch.matmul(delta_mat, GRN)
				# 保持delta_mat在gene_index列的值为init_col
				for i,idx in enumerate(gene_index):
					delta_mat[:, idx] = self.init_col[i].clone()

			# 模拟扰动后的mat，确保perturbed_mat的所有元素非负
			perturbed_mat = torch.clamp(self.wt_mat_t + delta_mat, min=0)
			
			# 趁着t转置之前顺便算了cos
			cos_avg_corr = cosine_similarity(perturbed_mat, ko_mat_t, 'mean').item()

			perturbed_mat = perturbed_mat.t()
			
			# 计算MSE损失
			mse_loss = nn.functional.mse_loss(perturbed_mat, ko_mat)
			
			## ps相关性、Spearman秩相关系数、余弦相似度、euclidean_distance、曼哈顿距离等
			pearson_avg_corr = 0
			spearman_avg_corr = 0 
			
			pearson_corr = pearson_corrcoef(perturbed_mat, ko_mat)
			spearman_corr = spearman_corrcoef(perturbed_mat, ko_mat)

			pearson_avg_corr = pearson_corr.mean().item()
			spearman_avg_corr = spearman_corr.mean().item()
			
			euclidean_avg_dist = 0
			manhattan_avg_dist = 0
			for i in range(perturbed_mat.shape[1]):
				perturb_col = perturbed_mat[:, i].unsqueeze(0)
				ko_col = ko_mat[:, i].unsqueeze(0)
				euclidean_avg_dist += pairwise_euclidean_distance(perturb_col,ko_col).item()
				manhattan_avg_dist += pairwise_manhattan_distance(perturb_col,ko_col).item()
			euclidean_avg_dist /= wt_mat.shape[1]
			manhattan_avg_dist /= wt_mat.shape[1]

			# 添加正则化项
			delta_GRN = GRN - GRN_init
			for i,idx in enumerate(gene_index):
				delta_GRN[idx] *= 2
				delta_GRN[:,idx] *= 2
			reg_loss = self.reg_weight * torch.norm(delta_GRN, p='fro') # F范数

			# 总损失函数
			total_loss = mse_loss  + reg_loss
			
			
			"""
			相关指标追加到对应的list
			"""
			# 记录损失
			mse_losses.append(mse_loss.item())
			reg_losses.append(reg_loss.item())
			total_losses.append(total_loss.item())
			
			# MAE
			mae_lists.append(mae_criterion(perturbed_mat, ko_mat).item())
			
			# RMSE
			RMSE_lists.append(math.sqrt(mse_loss.item()))
			
			# MAPE 
			MAPE_lists.append(mean_absolute_percentage_error(perturbed_mat, ko_mat).item())
			
			# 记录相关性指标
			pearson_avg_corrs.append(pearson_avg_corr)
			spearman_avg_corrs.append(spearman_avg_corr)
			cos_avg_corrs.append(cos_avg_corr)
			
			# 欧式距离
			euclidean_avg_dists.append(euclidean_avg_dist)
			
			# 曼哈顿距离
			manhattan_avg_dists.append(manhattan_avg_dist)
			
			# r2
			r2_score_list.append(r2_score(perturbed_mat, ko_mat).item())
			
			
			return total_loss
		
	# 各数据集不同，自行调整至较好效果
	reg_weight = 0.001 # 正则化权重参数
	lr = 0.001 # 学习率
	num_epochs = 100 # 设置训练轮数

	# 实例化损失函数
	loss_fn = MatrixRefinementLoss(init_col=init_col, wt_mat_t=wt_mat_t, reg_weight=reg_weight)

	# 定义优化器
	optimizer = optim.Adam([GRN], lr=lr)
	# optimizer = optim.RMSprop([GRN], lr=lr) # 另一个自适应学习率方法，有助于防止学习率过大而导致的震荡问题。
	# optimizer = optim.Adagrad([GRN], lr=lr) # 自适应学习率方法，对于稀疏梯度有良好的效果，但可能随着训练迭代导致学习率过早减小
	# optimizer = optim.Adamax([GRN], lr=lr) # Adam的一个变种，对非常稀疏的数据集有时会表现得更好。
	# optimizer = optim.Adadelta([GRN], lr=lr) # 另一种自适应学习率算法，不需要手动选择学习率，初始学习率为1比较好

	# 进行训练迭代直至收敛
	for epoch in range(num_epochs):
		optimizer.zero_grad()  # 清零梯度缓存

		total_epoch_loss = loss_fn()

		total_epoch_loss.backward()  # 反向传播计算梯度
		optimizer.step()  # 更新参数

		# 打印当前批次的损失值
		if not epoch % 50:
			print(
				f"{epoch} / {num_epochs} total loss:{total_losses[-1]:.4f} "
				f"mse: {mse_losses[-1]:.6f} "
				f"reg: {reg_losses[-1]:.5f} "
				f"ps: {pearson_avg_corrs[-1]:.4f} "
				f"spear: {spearman_avg_corrs[-1]:.4f} "
				f"cos: {cos_avg_corrs[-1]:.4f} "
				f"eucli: {euclidean_avg_dists[-1]:.4f} "
				f"mae: {mae_lists[-1]:.4f} "
				f"r2: {r2_score_list[-1]:.4f} "
				f"RMSE: {RMSE_lists[-1]:.4f} "
				f"MAPE: {MAPE_lists[-1]:.4f} "
				f"manhat: {manhattan_avg_dists[-1]:.4f} ")
    # 训练后的GRN矩阵
	GRN_refined = GRN.detach().clone()
	#可视化：
	# 增加全局字体大小
	font_size = 20  # 自定义一个较大的字体大小
	plt.rcParams.update({'font.size': font_size})
	# 绘制损失曲线
	plt.figure(figsize=(10, 6))

	plt.plot(range(1, num_epochs + 1), mse_losses, label='MSE Loss')
	plt.plot(range(1, num_epochs + 1), reg_losses, label='Regularization Loss')
	plt.plot(range(1, num_epochs + 1), total_losses, label='Total Loss')

	plt.xlabel('Epochs')
	plt.ylabel('Loss Value')
	plt.legend()
	plt.title('Losses During Training')
	plt.grid(True)
	plt.savefig(basePath + '/refine/Losses.png')
	# 绘制评价指标曲线
	plt.figure(figsize=(10, 6))

	plt.plot(range(1, num_epochs + 1), pearson_avg_corrs, label='Pearson corr')
	plt.plot(range(1, num_epochs + 1), spearman_avg_corrs, label='Spearman corr')
	plt.plot(range(1, num_epochs + 1), cos_avg_corrs, label='Cosine similarity')
	plt.plot(range(1, num_epochs + 1), r2_score_list, label='R2 score')

	plt.xlabel('Epochs')
	plt.ylabel('Value')
	plt.legend()
	plt.title('Performance Evaluation During Training')
	plt.grid(True)

	plt.savefig(basePath + '/refine/PerformanceEvaluation.png')
	fig, ax1 = plt.subplots(figsize=(12, 6))

	# 欧式距离曲线绘制于左侧Y轴
	color1 = '#0072BD'
	ax1.set_xlabel('Epochs')
	ax1.set_ylabel('Euclidean Avg Dists', color=color1)
	ax1.plot(range(1, num_epochs + 1), euclidean_avg_dists, color=color1)
	ax1.tick_params(axis='y', labelcolor=color1)

	# 曼哈顿距离曲线绘制于右侧Y轴
	ax2 = ax1.twinx()
	color2 = '#D95319'
	ax2.set_ylabel('Manhattan Avg Dists', color=color2)
	ax2.plot(range(1, num_epochs + 1), manhattan_avg_dists, color=color2)
	ax2.tick_params(axis='y', labelcolor=color2)

	# 共享X轴标签和标题
	plt.title('Distances Value During Training')
	plt.grid(True)

	# 不同Y轴可能需要手动调整刻度确保可读性
	# ax1.set_ylim(...)  # 调整欧式距离Y轴范围
	# ax2.set_ylim(...)  # 调整曼哈顿距离Y轴范围

	plt.legend(handles=[mpatches.Patch(color=color1, label='Euclidean Avg Dists'), 
						mpatches.Patch(color=color2, label='Manhattan Avg Dists')])
	plt.savefig(basePath + '/refine/Distances.png')
	# 将torch.Tensor转为numpy数组
	GRN_numpy = GRN_refined.cpu().numpy()

	# 将 lr 转换为科学计数法的字符串形式，不包括前导零
	lr_str = '{:.1e}'.format(lr)
	reg_str = '{:.1e}'.format(reg_weight)

	# 保存为.npy文件
	np.save(basePath + '/refine/refined_GRN.npy', GRN_numpy)
	
