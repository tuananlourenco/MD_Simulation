import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import argparse

# === Funções de equilibration.py (Mantidas como no código anterior) ===

def statisticalInefficiency_multiscale(A_n, B_n=None, fast=True, mintime=3):
    """
    Compute the (cross) statistical inefficiency of (two) timeseries using multiscale method from Chodera.
   
    Parameters
    ----------
    A_n : np.ndarray, float
        A_n[n] is nth value of timeseries A. Length is deduced from vector.
    B_n : np.ndarray, float, optional, default=None
        B_n[n] is nth value of timeseries B. Length is deduced from vector.
        If supplied, the cross-correlation of timeseries A and B will be estimated instead of the
        autocorrelation of timeseries A.
    fast : bool, optional, default=False
        f True, will use faster (but less accurate) method to estimate correlation
        time, described in Ref. [1] (default: False).
    mintime : int, optional, default=3
        minimum amount of correlation function to compute (default: 3)
        The algorithm terminates after computing the correlation time out to mintime when the
        correlation function first goes negative. Note that this time may need to be increased
        if there is a strong initial negative peak in the correlation function.

    Returns
    -------
    g : np.ndarray,
        g is the estimated statistical inefficiency (equal to 1 + 2 tau, where tau is the correlation time).
        We enforce g >= 1.0.

    Notes
    -----
    The same timeseries can be used for both A_n and B_n to get the autocorrelation statistical inefficiency.
    The fast method described in Ref [1] is used to compute g.

    References
    ----------
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
        histogram analysis method for the analysis of simulated and parallel tempering simulations.
        JCTC 3(1):26-41, 2007.

    Examples
    --------

    Compute statistical inefficiency of timeseries data with known correlation time.

    >>> from pymbar import testsystems
    >>> A_n = testsystems.correlated_timeseries_example(N=100000, tau=5.0)
    >>> g = statisticalInefficiency_multiscale(A_n, fast=True)

    """

    # Create numpy copies of input arguments.
    A_n = np.array(A_n)

    if B_n is not None:
        B_n = np.array(B_n)
    else:
        B_n = np.array(A_n)

    # Get the length of the timeseries.
    N = A_n.size

    # Be sure A_n and B_n have the same dimensions.
    if(A_n.shape != B_n.shape):
        raise Exception('A_n and B_n must have same dimensions.')

    # Initialize statistical inefficiency estimate with uncorrelated value.
    g = 1.0

    # Compute mean of each timeseries.
    mu_A = A_n.mean()
    mu_B = B_n.mean()

    # Make temporary copies of fluctuation from mean.
    dA_n = A_n.astype(np.float64) - mu_A
    dB_n = B_n.astype(np.float64) - mu_B

    # Compute estimator of covariance of (A,B) using estimator that will ensure C(0) = 1.
    sigma2_AB = (dA_n * dB_n).mean() # standard estimator to ensure C(0) = 1

    # Trap the case where this covariance is zero, and we cannot proceed.
    if(sigma2_AB == 0):
        # Original code raised 'ParameterError', which is not defined in this context.
        # Changed to a generic Exception for standalone execution.
        raise Exception('Sample covariance sigma_AB^2 = 0 -- cannot compute statistical inefficiency')

    # Accumulate the integrated correlation time by computing the normalized correlation time at
    # increasing values of t. Stop accumulating if the correlation function goes negative, since
    # this is unlikely to occur unless the correlation function has decayed to the point where it
    # is dominated by noise and indistinguishable from zero.
    t = 1
    increment = 1
    while (t < N - 1):

        # compute normalized fluctuation correlation function at time t
        C = np.sum(dA_n[0:(N - t)] * dB_n[t:N] + dB_n[0:(N - t)] * dA_n[t:N]) / (2.0 * float(N - t) * sigma2_AB)
        # Terminate if the correlation function has crossed zero and we've computed the correlation
        # function at least out to 'mintime'.
        if (C <= 0.0) and (t > mintime):
            break

        # Accumulate contribution to the statistical inefficiency.
        g += 2.0 * C * (1.0 - float(t) / float(N)) * float(increment)

        # Increment t and the amount by which we increment t.
        t += increment

        # Increase the interval if "fast mode" is on.
        if fast:
            increment += 1

    # g must be at least unity
    if (g < 1.0):
        g = 1.0

    # Return the computed statistical inefficiency.
    return g

def statisticalInefficiency_geyer(A_n, method='con'):
    """Compute the statistical inefficiency of a timeseries using the methods of Geyer.

    Parameters
    ----------
    A_n : np.ndarray, float
        A_n[n] is nth value of timeseries A. Length is deduced from vector.
    method : str, optional, default='con'
        The method to use; matches notation from `initseq` from `mcmc` R package by Geyer.
        'pos' : initial positive sequence (IPS) estimator
        'dec' : initial monotone sequence (IMS) estimator
        'con' : initial convex sequence (ICS) estimator

    Returns
    -------
    g : np.ndarray,
        g is the estimated statistical inefficiency (equal to 1 + 2 tau, where tau is the correlation time).
        We enforce g >= 1.0.

    Notes
    -----
    Implementation based on the `initseq` method from the `mcmc` R package by Geyer.

    References
    ----------
    [1] Geyer, CJ. Practical Markov chain Monte Carlo. Statistical Science 7(4):473-511, 1992.

    Examples
    --------

    Compute statistical inefficiency of timeseries data with known correlation time using different methods.

    >>> from pymbar import testsystems
    >>> A_n = testsystems.correlated_timeseries_example(N=100000, tau=10.0)
    >>> g_IPS = statisticalInefficiency_geyer(A_n, method='pos')
    >>> g_IMS = statisticalInefficiency_geyer(A_n, method='dec')
    >>> g_ICS = statisticalInefficiency_geyer(A_n, method='con')

    """

    if method not in ['pos', 'dec', 'con']:
        raise Exception("Unknown method '%s'; must be one of ['pos', 'dec', 'con']" % method)

    # Create numpy copies of input arguments.
    A_n = np.array(A_n)

    # Subtract off sample mean.
    A_n -= A_n.mean()

    # Get the length of the timeseries.
    N = A_n.size

    # Compute sample variance.
    gamma_zero = (A_n**2).sum() / N

    # Compute sequential covariance pairs.
    gamma_pos = list()
    for i in range(N//2): # Changed N/2 to N//2 for integer division
        lag1 = 2*i
        gam1 = (A_n[0:(N-lag1)] * A_n[lag1:N]).sum() / N
        lag2 = lag1 + 1
        gam2 = (A_n[0:(N-lag2)] * A_n[lag2:N]).sum() / N

        # Terminate if sum is no longer positive.
        if (gam1 + gam2) < 0.0:
            break

        # Otherwise, store the consecutive sum.
        gamma_pos.append(gam1+gam2)

    # Number of nonnegative values in array.
    ngamma = len(gamma_pos)

    # Compute IPS gamma sequence.
    gamma_pos = np.array(gamma_pos)

    # Compute IMS gamma sequence.
    gamma_dec = np.array(gamma_pos)
    for i in range(ngamma - 1):
        if gamma_dec[i] < gamma_dec[i+1]:
            gamma_dec[i] = gamma_dec[i+1]

    # Compute ICS gamma sequence.
    gamma_con = np.array(gamma_dec)
    for i in range(ngamma-1, 0, -1):
        gamma_con[i] -= gamma_con[i-1]

    # Pool adjacent violators (PAVA) algorithm.
    puff = np.zeros([ngamma], np.float64)
    nuff = np.zeros([ngamma], np.int32)
    nstep = 0
    for j in range(1,ngamma):
        puff[nstep] = gamma_con[j]
        nuff[nstep] = 1
        nstep += 1
        while (nstep > 1) and ((puff[nstep-1] / nuff[nstep-1]) < (puff[nstep-2] / nuff[nstep-2])):
            puff[nstep-2] += puff[nstep-1]
            nuff[nstep-2] += nuff[nstep-1]
            nstep -= 1

    j = 1
    for jstep in range(nstep):
        muff = puff[jstep] / nuff[jstep]
        for k in range(nuff[jstep]):
            gamma_con[j] = gamma_con[j-1] + muff
            j += 1

    # Compute sample variance estimates.
    var_pos = (2 * gamma_pos.sum() - gamma_zero) / N
    var_dec = (2 * gamma_dec.sum() - gamma_zero) / N
    var_con = (2 * gamma_con.sum() - gamma_zero) / N

    # Compute statistical inefficiencies from sample mean var = var(A_n) / (N/g)
    # g = var / (var(A_n)/N)
    var_uncorr = gamma_zero / N
    g_pos = var_pos / var_uncorr
    g_dec = var_dec / var_uncorr
    g_con = var_con / var_uncorr

    # Select appropriate g.
    if method == 'pos':
        g = g_pos
    elif method == 'dec':
        g = g_dec
    elif method == 'con':
        g = g_con

    # g must be at least unity
    if (g < 1.0):
        g = 1.0

    # Return the computed statistical inefficiency.
    return g

def detectEquilibration(A_t, nskip=1, method='multiscale'):
    """Automatically detect equilibrated region of a dataset using a heuristic that maximizes number of effectively uncorrelated samples.

    Parameters
    ----------
    A_t : np.ndarray
        timeseries
    nskip : int, optional, default=1
        number of samples to sparsify data by in order to speed equilibration detection
    method : str, optional, default='geyer'
        Method to use for computing statistical inefficiency; one of ['geyer', 'multiscale']

    Returns
    -------
    t : int
        start of equilibrated data
    g : float
        statistical inefficiency of equilibrated data
    Neff_max : float
        number of uncorrelated samples

    ToDo
    ----
    Consider implementing a binary search for Neff_max.

    Notes
    -----
    If your input consists of some period of equilibration followed by
    a constant sequence, this function treats the trailing constant sequence
    as having Neff = 1.

    Examples
    --------

    Determine start of equilibrated data for a correlated timeseries.

    >>> from pymbar import testsystems
    >>> A_t = testsystems.correlated_timeseries_example(N=1000, tau=10.0) # generate a test correlated timeseries
    >>> [t, g, Neff_max] = detectEquilibration(A_t) # compute indices of uncorrelated timeseries

    Determine start of equilibrated data for a correlated timeseries with a shift.

    >>> from pymbar import testsystems
    >>> A_t = testsystems.correlated_timeseries_example(N=1000, tau=10.0) + 2.0 # generate a test correlated timeseries
    >>> B_t = testsystems.correlated_timeseries_example(N=10000, tau=10.0) # generate a test correlated timeseries
    >>> C_t = np.concatenate([A_t, B_t])
    >>> [t, g, Neff_max] = detectEquilibration(C_t, nskip=50) # compute indices of uncorrelated timeseries

    """
    T = A_t.size

    # Special case if timeseries is constant.
    if A_t.std() == 0.0:
        return (0, 1, 1) # Changed from Neff=N to Neff=1 after issue #122

    if method == 'multiscale':
        statisticalInefficiency = statisticalInefficiency_multiscale
    elif method == 'geyer':
        statisticalInefficiency = statisticalInefficiency_geyer
    else:
        raise Exception("Method must be one of ['multiscale', 'geyer']")

    g_t = np.ones([T - 1], np.float32)
    Neff_t = np.ones([T - 1], np.float32)
    for t in range(0, T - 1, nskip):
        g_t[t] = statisticalInefficiency(A_t[t:T])
        Neff_t[t] = (T - t + 1) / g_t[t]
    Neff_max = Neff_t.max()
    t = Neff_t.argmax()
    g = g_t[t]

    return (t, g, Neff_max)


# === Código principal ===

# 1. Configurar os argumentos da linha de comando
parser = argparse.ArgumentParser(description="Analisa dados de simulação para detectar equilíbrio e plota os resultados.")
parser.add_argument("xvg_filename", help="Nome do arquivo .xvg de entrada (ex: energy.xvg)")
parser.add_argument("--line-color", default="tab:blue", help="Cor da linha da série temporal (padrão: tab:blue)")
parser.add_argument("--vline-color", default="red", help="Cor da linha de equilíbrio vertical (padrão: red)")
parser.add_argument("--xlim", nargs=2, type=float, help="Limites do eixo X (ex: --xlim 0 1000)")
parser.add_argument("--ylim", nargs=2, type=float, help="Limites do eixo Y (ex: --ylim 50 150)")
parser.add_argument("--legend-loc", default="best", help="Localização da legenda (ex: upper right, lower left, best)")
parser.add_argument("--xlabel", default="Tempo (ps)", help="Rótulo do eixo X")
parser.add_argument("--ylabel", default="Valor (ex: Densidade)", help="Rótulo do eixo Y")
parser.add_argument("--title-prefix", default="Detecção de Equilíbrio Automática", help="Prefixo do título do gráfico")
parser.add_argument("--full-equil-stats", "-f", action="store_true",
                    help="Calcular média e desvio padrão para *toda* a região equilibrada (padrão: apenas a segunda metade)")


args = parser.parse_args()

xvg_filename = args.xvg_filename
base_filename = os.path.splitext(xvg_filename)[0]
log_filename = f"{base_filename}-equilibrium-detect.log"

# Armazenar o stdout original para restaurar depois
original_stdout = sys.stdout

# Redirecionar a saída para o arquivo de log
try:
    with open(log_filename, 'w') as f:
        sys.stdout = f # Redireciona a saída padrão para o arquivo

        print(f"Iniciando análise para o arquivo: {xvg_filename}")

        # === 2. Carregue seu arquivo .xvg ===
        try:
            data = np.loadtxt(xvg_filename, comments=('@', '#'))
        except FileNotFoundError:
            print(f"Erro: O arquivo '{xvg_filename}' não foi encontrado. Por favor, certifique-se de que ele está no mesmo diretório que o script.")
            sys.exit(1)
        except Exception as e:
            print(f"Erro ao carregar o arquivo '{xvg_filename}': {e}")
            sys.exit(1)

        # === 3. Separe tempo e valores ===
        time = data[:, 0]
        values = data[:, 1]

        # === 4. Detecte a equilibração ===
        t_eq_idx, g, Neff = detectEquilibration(values)
        t_eq_time = time[t_eq_idx]

        print(f"Início da equilibração: índice {t_eq_idx} (tempo: {t_eq_time} ps)")
        print(f"Ineficiência estatística g: {g:.2f}")
        print(f"Número efetivo de amostras independentes: {Neff:.2f}")

        mean_equilibrated = None
        std_equilibrated = None
        stats_data_type = "" # Para usar no texto do log e do gráfico

        # === 5. Calcule a média e o desvio padrão dos dados equilibrados ===
        equilibrated_values = values[t_eq_idx:]
        if len(equilibrated_values) > 0:
            if args.full_equil_stats:
                data_for_stats = equilibrated_values
                stats_data_type = "Equilibrados (Total)"
            else:
                # Calcula o ponto de início da segunda metade
                start_second_half_idx = len(equilibrated_values) // 2
                data_for_stats = equilibrated_values[start_second_half_idx:]
                stats_data_type = "Equilibrados (2ª Metade)"

            if len(data_for_stats) > 0:
                mean_equilibrated = np.mean(data_for_stats)
                std_equilibrated = np.std(data_for_stats)
                print(f"Média dos dados {stats_data_type}: {mean_equilibrated:.4f}")
                print(f"Desvio padrão dos dados {stats_data_type}: {std_equilibrated:.4f}")
            else:
                print(f"Não há dados suficientes na região '{stats_data_type}' para calcular média e desvio padrão.")
        else:
            print("A região equilibrada não possui dados para calcular média e desvio padrão.")

        # === 6. Plotar resultado ===
        plt.figure(figsize=(10, 6))
        plt.plot(time, values, label="Série temporal", alpha=0.8, color=args.line_color)
        plt.axvline(t_eq_time, color=args.vline_color, linestyle='--', label=f'Equilíbrio: {t_eq_time:.1f} ps', linewidth=1.5)

        if args.xlim:
            plt.xlim(args.xlim[0], args.xlim[1])
        if args.ylim:
            plt.ylim(args.ylim[0], args.ylim[1])

        plt.xlabel(args.xlabel, fontsize=12)
        plt.ylabel(args.ylabel, fontsize=12)
        plt.legend(loc=args.legend_loc, fontsize=10)
        plt.title(f"{args.title_prefix} para {xvg_filename}", fontsize=14)
        plt.grid(True, linestyle=':', alpha=0.7)

        # Adicionar média e desvio padrão ao gráfico, se calculados
        if mean_equilibrated is not None and std_equilibrated is not None:
            xmin, xmax = plt.xlim()
            ymin, ymax = plt.ylim()

            text_x = xmax * 0.95
            text_y_mean = ymax - (ymax - ymin) * 0.05
            text_y_std = ymax - (ymax - ymin) * 0.12

            # Ajusta o rótulo no gráfico para ser mais conciso
            plot_stats_label = "Equil." if args.full_equil_stats else "Meia Equil."

            plt.text(text_x, text_y_mean, f'Média ({plot_stats_label}): {mean_equilibrated:.4f}',
                     fontsize=9, verticalalignment='top', horizontalalignment='right',
                     bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.5))
            plt.text(text_x, text_y_std, f'DP ({plot_stats_label}): {std_equilibrated:.4f}',
                     fontsize=9, verticalalignment='top', horizontalalignment='right',
                     bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.5))

        plt.tight_layout()
        plt.savefig(f"{base_filename}-equilibrium-plot.png", dpi=300)
        plt.show()

finally:
    sys.stdout = original_stdout
    print(f"Análise concluída. Resultados detalhados salvos em '{log_filename}'.")
    print(f"O gráfico foi salvo como '{base_filename}-equilibrium-plot.png' e exibido.")
