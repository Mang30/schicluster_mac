import time
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, diags, eye, save_npz
from scipy.sparse.linalg import norm
from scipy.ndimage import gaussian_filter
import cooler
import logging

# Import torch for MPS acceleration
try:
    import torch
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
    # Check if MPS is available on Mac
    MPS_AVAILABLE = torch.backends.mps.is_available()
    if MPS_AVAILABLE:
        logging.info("MPS acceleration is available and will be used for GPU acceleration.")
    else:
        logging.info("MPS acceleration is not available. Falling back to CPU.")
except ImportError:
    TORCH_AVAILABLE = False
    MPS_AVAILABLE = False
    logging.info("PyTorch not found. Using CPU-only implementation.")

# from ..cool import write_coo


def calc_sparsity(matrix):
    row, col = matrix.shape
    sparsity = matrix.nnz / row / col
    return sparsity


def scipy_to_torch_sparse(scipy_matrix, device):
    """Convert scipy sparse matrix to torch sparse tensor on specified device."""
    if not TORCH_AVAILABLE:
        raise RuntimeError("PyTorch is not available")
    
    coo = scipy_matrix.tocoo()
    indices = torch.from_numpy(np.vstack((coo.row, coo.col))).long()
    values = torch.from_numpy(coo.data).float()
    shape = coo.shape
    return torch.sparse_coo_tensor(indices, values, shape, device=device).coalesce()


def torch_sparse_to_scipy(torch_tensor):
    """Convert torch sparse tensor back to scipy sparse matrix."""
    torch_tensor = torch_tensor.coalesce().cpu()
    indices = torch_tensor.indices().numpy()
    values = torch_tensor.values().numpy()
    shape = torch_tensor.shape
    return csr_matrix((values, (indices[0], indices[1])), shape=shape)


def random_walk_mps(P, rp, tol):
    """MPS-accelerated random walk with restart using dense matrices."""
    if rp == 1:
        return P
    
    if not MPS_AVAILABLE:
        logging.warning("MPS not available, falling back to CPU")
        return random_walk_cpu(P, rp, tol)
    
    _start_time = time.time()
    device = torch.device("mps")
    n_genes = P.shape[0]
    
    # With 24GB GPU memory, we can handle very large matrices
    if n_genes > 50000:  # 大幅提高限制，支持超大矩阵
        logging.warning(f"Matrix size {n_genes}x{n_genes} is extremely large. Consider using CPU for memory safety.")
        return random_walk_cpu(P, rp, tol)
    
    try:
        # With 24GB GPU memory, we can be more aggressive
        logging.info(f"🚀 Using MPS for {n_genes}x{n_genes} matrix with 24GB GPU memory")
        
        # Convert sparse matrix to dense for MPS compatibility
        logging.info(f"Converting {n_genes}x{n_genes} sparse matrix to dense for MPS...")
        P_dense = torch.from_numpy(P.toarray()).float().to(device)
        I_dense = torch.eye(n_genes, device=device, dtype=torch.float32)
        
        # 计算内存使用
        matrix_memory_mb = (n_genes * n_genes * 4) / (1024 * 1024)  # float32 = 4 bytes
        total_memory_mb = matrix_memory_mb * 3  # P, I, Q 三个矩阵
        logging.info(f"💾 Estimated GPU memory usage: {total_memory_mb:.1f} MB")
        
        Q = P_dense.clone()
        
        for i in range(30):
            # Q_new = P @ (Q * (1 - rp) + rp * I)
            Q_scaled = Q * (1 - rp) + I_dense * rp
            Q_new = torch.mm(P_dense, Q_scaled)
            
            # Calculate convergence
            delta = torch.norm(Q - Q_new)
            
            # 释放旧的 Q 内存
            del Q
            Q = Q_new.clone()
            
            _end_time = time.time()
            logging.info(
                f'🔥 MPS Iter {i + 1} takes {(_end_time - _start_time):.3f} seconds. '
                f'Loss: {delta:.3f}')
            
            if delta < tol:
                logging.info(f"✅ Converged after {i + 1} iterations")
                break
        
        # Convert back to scipy sparse matrix
        logging.info("Converting MPS result back to sparse matrix...")
        Q_np = Q.cpu().numpy()
        
        # 清理 GPU 内存
        del P_dense, I_dense, Q
        if hasattr(torch.mps, 'empty_cache'):
            torch.mps.empty_cache()
            logging.info("🧹 GPU memory cleared")
        
        return csr_matrix(Q_np)
        
    except Exception as e:
        logging.warning(f"MPS acceleration failed: {e}. Falling back to CPU.")
        return random_walk_cpu(P, rp, tol)


def gaussian_filter_mps(matrix, std, pad):
    """MPS-accelerated Gaussian filter with 24GB memory optimization."""
    if not MPS_AVAILABLE:
        logging.warning("MPS not available for Gaussian filter, falling back to CPU")
        return gaussian_filter(matrix.astype(np.float32).toarray(), std, order=0, mode='mirror', truncate=pad)
    
    try:
        device = torch.device("mps")
        n_rows, n_cols = matrix.shape
        
        # 计算内存使用 (考虑输入矩阵、填充矩阵、卷积核、输出)
        matrix_memory_mb = (n_rows * n_cols * 4) / (1024 * 1024)  # float32
        kernel_size = int(2 * std * pad + 1)
        if kernel_size % 2 == 0:
            kernel_size += 1
        
        padding = kernel_size // 2
        padded_size = (n_rows + 2*padding) * (n_cols + 2*padding)
        total_memory_mb = (matrix_memory_mb + padded_size * 4 / (1024*1024)) * 2  # 输入+输出
        
        logging.info(f"🔥 MPS Gaussian filter: {n_rows}x{n_cols} matrix, 内存需求: {total_memory_mb:.1f}MB")
        
        # Convert to dense tensor on MPS
        dense_matrix = torch.from_numpy(matrix.toarray()).float().to(device)
        
        # Create Gaussian kernel
        # 1D Gaussian kernel
        x = torch.arange(kernel_size, device=device).float()
        x = x - (kernel_size - 1) / 2
        kernel_1d = torch.exp(-x**2 / (2 * std**2))
        kernel_1d = kernel_1d / kernel_1d.sum()
        
        # Create 2D kernel
        kernel_2d = kernel_1d.unsqueeze(0) * kernel_1d.unsqueeze(1)
        kernel_2d = kernel_2d.unsqueeze(0).unsqueeze(0)  # Add batch and channel dims
        
        # Apply convolution
        dense_matrix = dense_matrix.unsqueeze(0).unsqueeze(0)  # Add batch and channel dims
        
        # Pad with reflection
        dense_matrix = F.pad(dense_matrix, (padding, padding, padding, padding), mode='reflect')
        
        # Apply convolution
        filtered = F.conv2d(dense_matrix, kernel_2d, padding=0)
        filtered = filtered.squeeze(0).squeeze(0)  # Remove batch and channel dims
        
        # Return as numpy array
        result = filtered.cpu().numpy()
        
        # 清理 GPU 内存
        del dense_matrix, kernel_1d, kernel_2d, filtered
        if hasattr(torch.mps, 'empty_cache'):
            torch.mps.empty_cache()
        
        logging.info(f"✅ MPS Gaussian filter completed, GPU 内存已清理")
        return result
        
    except Exception as e:
        logging.warning(f"❌ MPS Gaussian filter failed: {e}. Falling back to CPU.")
        return gaussian_filter(matrix.astype(np.float32).toarray(), std, order=0, mode='mirror', truncate=pad)


def random_walk_cpu(P, rp, tol):
    """CPU version of random walk with restart."""
    if rp == 1:
        return P

    _start_time = time.time()
    n_genes = P.shape[0]
    I = eye(n_genes, dtype=np.float32)
    Q = P.copy()
    for i in range(30):
        Q_new = P.dot(Q * (1 - rp) + rp * I)
        delta = norm(Q - Q_new)
        Q = Q_new.copy()
        sparsity = calc_sparsity(Q)
        _end_time = time.time()
        logging.debug(
            f'CPU Iter {i + 1} takes {(_end_time - _start_time):.3f} seconds. '
            f'Loss: {delta:.3f}; Sparsity: {sparsity:.3f}')
        if delta < tol:
            break
    return Q


def random_walk_auto(P, rp, tol, use_mps=True, max_memory_gb=24):
    """
    智能选择 MPS 或 CPU，专为 24GB GPU 内存优化
    
    参数:
    - P: 转移概率矩阵 (scipy.sparse)
    - rp: 重启概率
    - tol: 收敛容差
    - use_mps: 是否启用 MPS 加速
    - max_memory_gb: 最大 GPU 内存 (默认 24GB)
    """
    n_genes = P.shape[0]
    sparsity = P.nnz / (n_genes * n_genes)
    
    # 计算内存需求 (P, I, Q 三个密集矩阵)
    memory_needed_gb = (n_genes * n_genes * 4 * 3) / (1024**3)
    
    # 24GB GPU 可以处理非常大的矩阵，极低的阈值
    mps_threshold = 100  # 极低阈值，几乎所有矩阵都用 MPS
    
    # 决策逻辑
    use_mps_decision = (
        use_mps and 
        MPS_AVAILABLE and 
        n_genes >= mps_threshold and
        memory_needed_gb <= max_memory_gb and
        sparsity > 0.0001  # 极稀疏矩阵仍用 CPU
    )
    
    if use_mps_decision:
        logging.info(f"🚀 MPS 加速: {n_genes}x{n_genes} 矩阵 (内存需求: {memory_needed_gb:.2f}GB, 稀疏度: {sparsity:.4f})")
        return random_walk_mps(P, rp, tol)
    else:
        if use_mps and MPS_AVAILABLE:
            if memory_needed_gb > max_memory_gb:
                reason = f"内存超限 ({memory_needed_gb:.2f}GB > {max_memory_gb}GB)"
            elif n_genes < mps_threshold:
                reason = f"矩阵太小 ({n_genes} < {mps_threshold})"
            elif sparsity <= 0.0001:
                reason = f"过于稀疏 (稀疏度: {sparsity:.6f})"
            else:
                reason = "未知原因"
            logging.info(f"💻 CPU 回退: {n_genes}x{n_genes} 矩阵 ({reason})")
        else:
            logging.info(f"💻 CPU 处理: {n_genes}x{n_genes} 矩阵 (MPS 未启用或不可用)")
        return random_walk_cpu(P, rp, tol)


def impute_chromosome(chrom,
                      resolution,
                      output_path,
                      scool_url=None,
                      contact_path=None,
                      chrom_size_path=None,
                      logscale=False,
                      pad=1,
                      std=1,
                      rp=0.5,
                      tol=0.01,
                      window_size=5000000000,
                      step_size=10000000,
                      output_dist=5000000000,
                      min_cutoff=0,
                      chrom1=1,
                      pos1=2,
                      chrom2=5,
                      pos2=6,
                      use_mps=True):
    """
    Impute chromosome contact matrix with optional MPS GPU acceleration.
    
    Parameters:
    -----------
    use_mps : bool, default=True
        Whether to use MPS (Metal Performance Shaders) acceleration on Mac.
        If False or MPS not available, falls back to CPU computation.
    """
    # Log acceleration status
    if use_mps:
        if MPS_AVAILABLE:
            logging.info(f"Starting chromosome {chrom} imputation with MPS GPU acceleration")
        else:
            logging.info(f"MPS requested but not available. Using CPU for chromosome {chrom}")
    else:
        logging.info(f"Starting chromosome {chrom} imputation with CPU (MPS disabled)")
    
    if scool_url is not None:
        cell_cool = cooler.Cooler(scool_url)
        A = cell_cool.matrix(balance=False, sparse=True).fetch(chrom)
        # A = A + diags(A.diagonal())
        n_bins = A.shape[0]
    elif contact_path is not None:
        if chrom_size_path is not None:
            chrom_sizes = pd.read_csv(chrom_size_path, sep='\t', index_col=0, header=None).squeeze(axis=1)
            if chrom=='all':
                from schicluster.cool.utilities import get_chrom_offsets
                bins_df = cooler.binnify(chrom_sizes, resolution)
                chrom_offset = get_chrom_offsets(bins_df)
                n_bins = bins_df.shape[0]
            else:
                n_bins = (chrom_sizes.loc[chrom] // resolution) + 1
        else:
            print("ERROR : Must provide chrom_size_path if using contact file as input")
            return
        A = pd.read_csv(contact_path, sep='\t', header=None, index_col=None, comment='#')[[chrom1, pos1, chrom2, pos2]]
        if chrom=='all':
            A = A[A[chrom1].isin(chrom_offset) & A[chrom2].isin(chrom_offset)]
            A[pos1] = A[chrom1].map(chrom_offset) + (A[pos1] - 1) // resolution
            A[pos2] = A[chrom2].map(chrom_offset) + (A[pos2] - 1) // resolution
        else:
            A = A.loc[(A[chrom1]==chrom) & (A[chrom2]==chrom)]
            A[[pos1, pos2]] = (A[[pos1, pos2]] - 1) // resolution
        A = A.groupby(by=[pos1, pos2])[chrom1].count().reset_index()
        A = csr_matrix((A[chrom1].astype(np.int32), (A[pos1], A[pos2])), (n_bins, n_bins))
        A = A + A.T
    else:
        print("ERROR : Must provide either scool_url or contact_file_path")
        return

    ws = int(window_size // resolution)
    ss = int(step_size // resolution)

    # log transform
    if logscale:
        A.data = np.log2(A.data + 1)

    # Remove diagonal before convolution
    A = A - diags(A.diagonal())

    # Gaussian convolution
    start_time = time.time()
    if pad > 0:
        # Use MPS accelerated Gaussian filter if available and requested
        if use_mps and MPS_AVAILABLE:
            logging.info("Using MPS-accelerated Gaussian convolution")
            A_filtered = gaussian_filter_mps(A, std, pad)
            A = csr_matrix(A_filtered)
        else:
            logging.info("Using CPU Gaussian convolution")
            # full matrix step
            A = gaussian_filter(A.astype(np.float32).toarray(),
                                std, order=0, mode='mirror', truncate=pad)
            A = csr_matrix(A)
    # else:
    #     A = A + A.T
    end_time = time.time()
    logging.debug(f'Convolution takes {end_time - start_time:.3f} seconds')

    # Remove diagonal before RWR
    A = A - diags(A.diagonal())

    # Random Walk with Restart
    start_time = time.time()
    acceleration_method = "MPS" if (use_mps and MPS_AVAILABLE) else "CPU"
    logging.info(f"Using {acceleration_method} acceleration for Random Walk with Restart")
    
    if ws >= n_bins or rp == 1:
        B = A + diags((A.sum(axis=0).A.ravel() == 0).astype(int))
        d = diags(1 / B.sum(axis=0).A.ravel())
        P = d.dot(B).astype(np.float32)
        E = random_walk_auto(P, rp, tol, use_mps)
    else:
        # if the chromosome is too large, compute by chunks
        idx = (np.repeat(np.arange(ws), ws), np.tile(np.arange(ws), ws))
        idxfilter = (np.abs(idx[1] - idx[0]) < (output_dist // resolution + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
        # first filter
        idxfilter = ((idx[0] + idx[1]) < (ws + ss))
        idx1 = (idx[0][idxfilter], idx[1][idxfilter])
        mask1 = csr_matrix((np.ones(len(idx1[0])), (idx1[0], idx1[1])),
                           (ws, ws))
        # last filter
        idxfilter = ((idx[0] + idx[1]) >= (
                (n_bins - ws) // ss * 2 + 1) * ss + 3 * ws - 2 * n_bins)
        idx2 = (idx[0][idxfilter], idx[1][idxfilter])
        mask2 = csr_matrix((np.ones(len(idx2[0])), (idx2[0], idx2[1])),
                           (ws, ws))
        # center filter
        idxfilter = np.logical_and((idx[0] + idx[1]) < (ws + ss),
                                   (idx[0] + idx[1]) >= (ws - ss))
        idx0 = (idx[0][idxfilter], idx[1][idxfilter])
        mask0 = csr_matrix((np.ones(len(idx0[0])), (idx0[0], idx0[1])),
                           (ws, ws))

        start_time = time.time()
        E = csr_matrix(A.shape, dtype=np.float32)
        for ll in [x for x in range(0, n_bins - ws, ss)] + [n_bins - ws]:
            B = A[ll:(ll + ws), ll:(ll + ws)]
            B = B + diags((B.sum(axis=0).A.ravel() == 0).astype(int))
            d = diags(1 / B.sum(axis=0).A.ravel())
            P = d.dot(B).astype(np.float32)
            Etmp = random_walk_auto(P, rp, tol, use_mps)
            if ll == 0:
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask1)
            elif ll == (n_bins - ws):
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask2)
            else:
                E[ll:(ll + ws), ll:(ll + ws)] += Etmp.multiply(mask0)
    logging.debug(f'RWR takes {time.time() - start_time:.3f} seconds')

    # Normalize
    start_time = time.time()
    E += E.T
    d = E.sum(axis=0).A.ravel()
    d[d == 0] = 1
    b = diags(1 / np.sqrt(d))
    E = b.dot(E).dot(b)
    logging.debug(f'SQRTVC takes {time.time() - start_time:.3f} seconds')

    start_time = time.time()
    # mask the lower triangle of E
    # TODO This part is MEM intensive, the mask below can be combined with the chunk mask above
    idx = np.triu_indices(E.shape[0], 0)
    if (output_dist // resolution + 1) < n_bins:
        # longest distance filter mask
        idxfilter = ((idx[1] - idx[0]) < (output_dist // resolution + 1))
        idx = (idx[0][idxfilter], idx[1][idxfilter])
    mask = csr_matrix((np.ones(len(idx[0])), (idx[0], idx[1])),
                      E.shape,
                      dtype=np.float32)
    E = E.tocsr().multiply(mask)
    logging.debug(f'Filter takes {time.time() - start_time:.3f} seconds')

    # TODO put this part inside RWR, before normalize
    # min_cutoff = tol/
    # Make values < min_cutoff to 0
    if min_cutoff > 0:
        s_before = calc_sparsity(E)
        E = E.multiply(E > min_cutoff)
        s_after = calc_sparsity(E)
        logging.debug(f'Mask values smaller than {min_cutoff}. Sparsity before {s_before:.3f}, after {s_after:.3f}')

    # save to file
    # write_coo(output_path, E)
    save_npz(output_path, E)

    return
