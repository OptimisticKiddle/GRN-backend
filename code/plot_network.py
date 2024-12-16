import scipy.io
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import pandas as pd
import colorsys
import random
import json
import scipy.io
import numpy as np
import networkx as nx
import plotly.graph_objects as go
import pandas as pd
import colorsys
import random
import json

def generate_pastel_colors(n):
    """生成柔和的彩色色板"""
    # 设置固定的随机种子
    random.seed(42)
    
    colors = []
    for i in range(n):
        hue = i/n
        saturation = 0.3
        value = 0.95
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        color = '#{:02x}{:02x}{:02x}'.format(
            int(rgb[0]*255), 
            int(rgb[1]*255), 
            int(rgb[2]*255)
        )
        colors.append(color)
    # 移除随机打乱
    # random.shuffle(colors)
    return colors

def load_grn_data(mm_file, csv_file):
    adj_matrix = scipy.io.mmread(mm_file)
    adj_matrix = adj_matrix.tocoo()
    genes = pd.read_csv(csv_file, header=0).iloc[:,0].tolist()
    return adj_matrix, genes

def create_network(adj_matrix, genes, threshold=0.5):
    G = nx.Graph()
    
    # 添加节点
    for gene in genes:
        G.add_node(gene)
    
    # 添加边
    for i, j, data in zip(adj_matrix.row, adj_matrix.col, adj_matrix.data):
        if i < j:  # 避免重复边
            weight = abs(data)
            if weight >= threshold:
                G.add_edge(genes[i], genes[j], weight=weight)
    
    return G

def save_node_properties(pos, node_sizes, node_colors, file_path):
    """保存节点属性到JSON文件"""
    node_properties = {
        'positions': {node: {'x': float(coord[0]), 'y': float(coord[1])} 
                     for node, coord in pos.items()},
        'sizes': {node: float(size) for node, size in node_sizes.items()},
        'colors': node_colors
    }
    
    with open(file_path, 'w') as f:
        json.dump(node_properties, f)

def load_node_properties(file_path):
    """从JSON文件加载节点属性"""
    with open(file_path, 'r') as f:
        properties = json.load(f)
    
    # 将位置信息转换回numpy数组格式
    positions = {node: np.array([coord['x'], coord['y']]) 
                for node, coord in properties['positions'].items()}
    
    return positions, properties['sizes'], properties['colors']

def plot_network(G, min_edge_weight=0.5, max_nodes=100, node_properties=None):
    # 设置固定的随机种子
    random.seed(42)
    np.random.seed(42)
    
    if len(G.nodes()) > max_nodes:
        edges = sorted(G.edges(data=True), key=lambda x: x[2]['weight'], reverse=True)
        important_nodes = set()
        important_edges = []
        
        for edge in edges:
            if len(important_nodes) >= max_nodes:
                break
            important_nodes.add(edge[0])
            important_nodes.add(edge[1])
            important_edges.append(edge)
        
        G = nx.Graph()
        G.add_nodes_from(important_nodes)
        G.add_weighted_edges_from([(e[0], e[1], e[2]['weight']) for e in important_edges])
    
    if node_properties is None:
        # 如果没有提供节点属性，则重新计算
        pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), iterations=50, seed=42)
        colors = generate_pastel_colors(len(G.nodes()))
        node_colors = dict(zip(G.nodes(), colors))
        
        degrees = dict(G.degree())
        max_degree = max(degrees.values())
        node_sizes = {node: 20 + (degrees[node] / max_degree) * 80 for node in G.nodes()}
    else:
        # 使用提供的节点属性
        pos, node_sizes, node_colors = node_properties
        
        # 对于新出现的节点，计算新的属性
        new_nodes = set(G.nodes()) - set(pos.keys())
        if new_nodes:
            new_pos = nx.spring_layout(G.subgraph(new_nodes), k=1/np.sqrt(len(G.nodes())), iterations=50, seed=42)
            pos.update(new_pos)
            
            degrees = dict(G.degree())
            max_degree = max(degrees.values())
            new_sizes = {node: 20 + (degrees[node] / max_degree) * 80 for node in new_nodes}
            node_sizes.update(new_sizes)
            
            new_colors = generate_pastel_colors(len(new_nodes))
            new_node_colors = dict(zip(new_nodes, new_colors))
            node_colors.update(new_node_colors)
    
    # 创建节点轨迹
    node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers+text',
        textposition='top center',
        hoverinfo='text',
        marker=dict(
            size=[],
            color=[],
            line=dict(color='white', width=0.5)
        )
    )
    
    # 添加节点
    for node in G.nodes():
        x, y = pos[node]
        node_trace['x'] += tuple([x])
        node_trace['y'] += tuple([y])
        node_trace['text'] += tuple([node])
        node_trace['marker']['size'] += tuple([node_sizes[node]])
        node_trace['marker']['color'] += tuple([node_colors[node]])
    
    # 创建边轨迹
    edge_traces = []
    for edge in G.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        weight = edge[2]['weight']
        
        edge_trace = go.Scatter(
            x=[x0, x1],
            y=[y0, y1],
            mode='lines',
            line=dict(
                width=weight * 2,
                color='rgba(169,169,169,0.5)'
            ),
            hoverinfo='text',
            text=f'{edge[0]} - {edge[1]}<br>Weight: {weight:.3f}'
        )
        edge_traces.append(edge_trace)
    
    # 创建图形
    fig = go.Figure(data=edge_traces + [node_trace],
                   layout=go.Layout(
                       title='Gene Regulatory Network',
                       titlefont=dict(size=16),
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20,l=5,r=5,t=40),
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       plot_bgcolor='white',
                       width=1000,
                       height=1000
                   ))
    
    return fig, (pos, node_sizes, node_colors)

def visualize_grn(mm_file, csv_file, threshold=0.5, max_nodes=100, node_properties=None, save_properties=False):
    try:
        adj_matrix, genes = load_grn_data(mm_file, csv_file)
        print("Successfully loaded data")
        print(f"Matrix shape: {adj_matrix.shape}")
        print(f"Number of genes: {len(genes)}")
        
        G = create_network(adj_matrix, genes, threshold)
        print(f"Created network with {len(G.nodes())} nodes and {len(G.edges())} edges")
        
        fig, properties = plot_network(G, threshold, max_nodes, node_properties)
        print("Network visualization completed")
        
        if save_properties:
            properties_file = f"../static/{GSE_id}/{GSM_id}/refine/node_properties.json"
            save_node_properties(*properties, properties_file)
            print(f"Node properties saved to: {properties_file}")
        
        return fig, properties
        
    except Exception as e:
        print(f"Error during visualization: {str(e)}")
        return None, None

    
def create_network_from_numpy(grn_matrix, genes, threshold=0.5):
    G = nx.Graph()
    
    # 添加节点
    for gene in genes:
        G.add_node(gene)
    
    # 规范化权重
    max_weight = np.max(np.abs(grn_matrix))
    normalized_matrix = np.abs(grn_matrix) / max_weight
    
    # 只处理上三角矩阵
    n = len(genes)
    for i in range(n):
        for j in range(i+1, n):
            weight = normalized_matrix[i, j]
            if weight >= threshold:
                G.add_edge(genes[i], genes[j], weight=weight)
    
    return G

def plot_refined_network(G, max_nodes=100, lr_str='default', reg_str='default'):
    if len(G.nodes()) > max_nodes:
        edges = sorted(G.edges(data=True), key=lambda x: x[2]['weight'], reverse=True)
        important_nodes = set()
        important_edges = []
        
        for edge in edges:
            if len(important_nodes) >= max_nodes:
                break
            important_nodes.add(edge[0])
            important_nodes.add(edge[1])
            important_edges.append(edge)
        
        G = nx.Graph()
        G.add_nodes_from(important_nodes)
        G.add_weighted_edges_from([(e[0], e[1], e[2]['weight']) for e in important_edges])
    
    # 读取保存的节点属性
    properties_file = f"../static/{GSE_id}/{GSM_id}/refine/node_properties.json"
    try:
        with open(properties_file, 'r') as f:
            saved_properties = json.load(f)
        pos = {node: [coord[0], coord[1]] for node, coord in saved_properties['positions'].items()}
        node_colors = saved_properties['colors']
        node_sizes = saved_properties['sizes']
    except:
        # 如果无法读取保存的属性，使用默认的布局
        pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), iterations=50)
        colors = generate_pastel_colors(len(G.nodes()))
        node_colors = dict(zip(G.nodes(), colors))
        degrees = dict(G.degree())
        max_degree = max(degrees.values()) if degrees else 1
        node_sizes = {node: 20 + (degrees[node] / max_degree) * 80 for node in G.nodes()}
    
    node_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='markers+text',
        textposition='top center',
        hoverinfo='text',
        marker=dict(
            size=[],
            color=[],
            line=dict(color='white', width=0.5)
        )
    )
    
    for node in G.nodes():
        x, y = pos[node]
        node_trace['x'] += tuple([x])
        node_trace['y'] += tuple([y])
        node_trace['text'] += tuple([node])
        node_trace['marker']['size'] += tuple([node_sizes[node]])
        node_trace['marker']['color'] += tuple([node_colors[node]])
    
    edge_traces = []
    for edge in G.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        weight = edge[2]['weight']
        
        edge_trace = go.Scatter(
            x=[x0, x1],
            y=[y0, y1],
            mode='lines',
            line=dict(
                width=weight * 2,
                color='rgba(169,169,169,0.5)'
            ),
            hoverinfo='text',
            text=f'{edge[0]} - {edge[1]}<br>Weight: {weight:.3f}'
        )
        edge_traces.append(edge_trace)
    
    fig = go.Figure(data=edge_traces + [node_trace],
                   layout=go.Layout(
                       title=f'Refined Gene Regulatory Network (lr={lr_str}, reg={reg_str})',
                       titlefont=dict(size=16),
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20,l=5,r=5,t=40),
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       plot_bgcolor='white',
                       width=1000,
                       height=1000
                   ))
    
    return fig

def visualize_network(GRN_matrix, genes, output_path, lr_str='default', reg_str='default'):
    G = create_network_from_numpy(GRN_matrix, genes, threshold=0.5)
    fig = plot_refined_network(G, max_nodes=100, lr_str=lr_str, reg_str=reg_str)
    
    output_html = f"{output_path}_network.html"
    fig.write_html(output_html)

def get_properties(gse_id:str,gsm_id:str):
    global GSE_id
    global GSM_id
    GSE_id = 'GSE' + gse_id
    GSM_id = 'GSM' + gsm_id
    base_path = f'../static/{GSE_id}/{GSM_id}'
     # 生成baseGRN并保存节点属性
    base_mm_file = f"{base_path}/GRN.mm"
    base_csv_file = f"{base_path}/GRN_row_names.csv"
    
    visualize_grn(
        mm_file=base_mm_file,
        csv_file=base_csv_file,
        threshold=0.5,
        max_nodes=100,
        save_properties=True
    )

def get_refined_network(gse_id:str,gsm_id:str):
    GSE_id = 'GSE' + gse_id
    GSM_id = 'GSM' + gsm_id
    base_path = f'../static/{GSE_id}/{GSM_id}'
    learning_rate = 0.001
    reg_lambda = 0.01
    # 读取数据文件
    wt_csv_path = f"{base_path}/WT/wt_avg_expr_by_cluster.csv"
    ko_csv_path = f"{base_path}/KO/ko_avg_expr_by_cluster.csv"

    wt_data = pd.read_csv(wt_csv_path, index_col=0)
    ko_data = pd.read_csv(ko_csv_path, index_col=0)

    genes = wt_data.index.tolist()
    output_base_path = f"{base_path}/refine"

    GRN_path = f"{base_path}/refine/refined_GRN.npy"
    GRN = np.load(GRN_path)

    visualize_network(
		GRN_matrix=GRN,
		genes=genes,
		output_path=f"{output_base_path}/refined_GRN",
		lr_str=str(learning_rate),
		reg_str=str(reg_lambda)
    )

    