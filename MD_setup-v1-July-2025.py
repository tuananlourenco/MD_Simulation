import subprocess
import configparser
import os
import re
import shutil # Importado para operações de cópia de arquivos (shutil.copy2)

def clone_repository(repo_url, clone_dir="MD_Simulation"):
    """
    Clona um repositório Git para um diretório especificado.
    Se o diretório de destino já existir e não estiver vazio, a clonagem é pulada.

    Args:
        repo_url (str): A URL do repositório Git a ser clonado.
        clone_dir (str): O nome do diretório onde o repositório será clonado.
                         Por padrão, cria uma pasta chamada "MD_Simulation" no diretório atual.
    Returns:
        bool: True se o repositório foi clonado ou já existia, False em caso de erro.
    """
    # NOVO: Altera o URL de SSH para HTTPS
    # Se o URL for do tipo 'git@github.com:user/repo.git', converte para 'https://github.com/user/repo.git'
    if repo_url.startswith("git@github.com:"):
        repo_url = repo_url.replace("git@github.com:", "https://github.com/").replace(".git", ".git") # Mantém .git
        print(f"URL de clonagem ajustado para HTTPS: '{repo_url}'")


    target_path = os.path.join(os.getcwd(), clone_dir)

    if os.path.exists(target_path) and os.listdir(target_path):
        print(f"O repositório '{repo_url}' parece já ter sido clonado em '{target_path}'. Pulando a clonagem.")
        return True

    print(f"Clonando o repositório '{repo_url}' para '{target_path}'...")
    try:
        # Adicione o argumento '--depth 1' para clonar apenas o histórico mais recente, se o repositório for grande.
        subprocess.run(["git", "clone", repo_url, clone_dir], check=True, stderr=subprocess.PIPE, text=True)
        print("Repositório clonado com sucesso!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Erro ao clonar o repositório: {e.returncode}")
        if e.stderr:
            print(f"Detalhes do erro: {e.stderr}")
        else:
            print("Nenhum detalhe de erro disponível de stderr.")
        return False
    except FileNotFoundError:
        print("Erro: O comando 'git' não foi encontrado.")
        print("Certifique-se de que o Git está instalado e acessível no PATH do sistema.")
        return False
    except Exception as e:
        print(f"Um erro inesperado ocorreu durante a clonagem: {e}")
        return False

def read_user_variables(file_path):
    """
    Lê o arquivo user_variables.txt e armazena as variáveis em um objeto ConfigParser.
    Ignora comentários (linhas começando com # ou ;) e linhas vazias.

    Args:
        file_path (str): O caminho para o arquivo user_variables.txt.

    Returns:
        configparser.ConfigParser: Um objeto ConfigParser contendo as variáveis.
                                   Retorna None se o arquivo não for encontrado ou houver erro.
    """
    user_vars = configparser.ConfigParser(allow_no_value=True)
    try:
        read_files = user_vars.read(file_path)
        if not read_files:
            print(f"Aviso: O arquivo '{file_path}' foi encontrado, mas nenhuma seção válida foi lida.")
            print("Verifique se o arquivo não está vazio ou contém apenas comentários/erros de formatação.")
            return None

        print(f"Variáveis lidas de '{file_path}' com sucesso!")
        print(f"Seções encontradas pelo ConfigParser: {user_vars.sections()}")

        # ConfigParser converte nomes de seções para minúsculas. Verificamos 'general'.
        if 'general' not in user_vars.sections():
            print("\nERRO DE DEPURACÃO: A seção '[General]' (ou 'general' em minúsculas) NÃO foi encontrada pelo ConfigParser.")
            print("Por favor, verifique a formatação do seu user_variables.txt, especialmente a linha '[General]'.")
            print("Pode haver caracteres invisíveis ou formatação incorreta (ex: espaços extras antes/depois dos colchetes ou do nome da seção).")
            return None

        return user_vars
    except configparser.Error as e:
        print(f"Erro ao parsear o arquivo de variáveis: {e}")
        return None
    except FileNotFoundError:
        print(f"Erro: O arquivo '{file_path}' não foi encontrado.")
        return None

def create_main_simulation_directory(config_data):
    """
    Cria o diretório principal para as simulações com base nas variáveis do user_variables.
    O nome do diretório será Project_{overleaf_name}_{start_year}_MD_simulations.

    Args:
        config_data (configparser.ConfigParser): Objeto ConfigParser com as variáveis do usuário.

    Returns:
        str: O caminho completo para o diretório criado, ou None em caso de erro.
    """
    try:
        # Acessando as variáveis da seção 'General' para o nome do diretório principal
        # Usamos .get() com um valor padrão para evitar KeyError se a variável estiver faltando.
        overleaf_name = config_data['general'].get('overleaf_name', 'UnnamedProject')
        start_year = config_data['general'].get('start_year', 'YYYY')

        # Montando o nome do diretório principal
        dir_name = f"Project_{overleaf_name}_{start_year}_MD_simulations"
        
        # O diretório será criado no mesmo nível do script
        full_path = dir_name 
        
        os.makedirs(full_path, exist_ok=True)
        print(f"\nDiretório principal de simulação criado: '{full_path}'")
        return full_path

    except KeyError as e:
        print(f"Erro: Variável '{e}' necessária para o nome do diretório principal não encontrada ou vazia.")
        print("Por favor, verifique se 'overleaf_name' e 'start_year' estão definidos E PREENCHIDOS na seção [General] de user_variables.txt.")
        return None
    except Exception as e:
        print(f"Um erro inesperado ocorreu ao criar o diretório principal: {e}")
        return None

def create_run_subdirectories_actual(run_parent_dir, config_data):
    """
    Cria os 5 subdiretórios específicos de etapas da simulação dentro do diretório RUN.
    Retorna o caminho do diretório 'min_steep_10_500000' se criado com sucesso, senão None.

    Args:
        run_parent_dir (str): O caminho para o diretório RUN.
        config_data (configparser.ConfigParser): Objeto ConfigParser com as variáveis do usuário.
    Returns:
        str: O caminho completo para o diretório de minimização se bem-sucedido, senão None.
    """
    if not run_parent_dir:
        print("Não foi possível criar diretórios de etapas: o diretório RUN não foi especificado.")
        return None # Indica falha

    print(f"\nCriando diretórios de etapas da simulação em '{run_parent_dir}'...")

    min_dir_path = None # Variável para armazenar o caminho do diretório de minimização

    try:
        # 1. Energy Minimization
        integrator_em = config_data['energy minimization']['integrator']
        emtol_em = config_data['energy minimization']['emtol']
        nsteps_em = config_data['energy minimization']['nsteps']
        em_dir_name = f"min_{integrator_em}_{emtol_em}_{nsteps_em}"
        em_path = os.path.join(run_parent_dir, em_dir_name)
        os.makedirs(em_path, exist_ok=True)
        print(f"  Criado: '{em_path}'")
        min_dir_path = em_path # Armazena o caminho do diretório de minimização

        # 2. NVT 2 stage
        dt_nvt2 = float(config_data['nvt 2 stage']['dt'])
        nsteps_nvt2 = int(config_data['nvt 2 stage']['nsteps'])
        time_nvt2 = (dt_nvt2 * nsteps_nvt2) / 1000 # Tempo em ns
        time_nvt2_formatted = f"{time_nvt2:.2f}".replace('.', 'p') # Substitui '.' por 'p' para nomes de arquivos
        t_initial_nvt2 = config_data['nvt 2 stage']['t-initial']
        tcoupl_nvt2 = config_data['nvt 2 stage']['tcoupl']
        nvt2_dir_name = f"NVT_thermalization_{dt_nvt2}_{time_nvt2_formatted}ns_{nsteps_nvt2}_{t_initial_nvt2}-K_{tcoupl_nvt2}"
        nvt2_path = os.path.join(run_parent_dir, nvt2_dir_name)
        os.makedirs(nvt2_path, exist_ok=True)
        print(f"  Criado: '{nvt2_path}'")

        # 3. NPT Stage 3 (Annealing) - ANTERIORMENTE NPT Stage 2
        dt_npt3 = float(config_data['npt stage 3']['dt']) # Alterado de npt2 para npt3
        nsteps_npt3 = int(config_data['npt stage 3']['nsteps']) # Alterado de npt2 para npt3
        time_npt3 = (dt_npt3 * nsteps_npt3) / 1000 # Tempo em ns # Alterado de npt2 para npt3
        time_npt3_formatted = f"{time_npt3:.2f}".replace('.', 'p') # Alterado de npt2 para npt3
        t_initial_npt3 = config_data['npt stage 3']['t-initial'] # Alterado de npt2 para npt3
        t_final_npt3 = config_data['npt stage 3']['t-final'] # Alterado de npt2 para npt3
        tcoupl_npt3 = config_data['npt stage 3']['tcoupl'] # Alterado de npt2 para npt3
        ref_p_npt3 = config_data['npt stage 3']['ref-p'] # Alterado de npt2 para npt3
        pcoupl_npt3 = config_data['npt stage 3']['pcoupl'] # Alterado de npt2 para npt3
        npt3_dir_name = f"NPT_annealing_{dt_npt3}_{time_npt3_formatted}ns_{nsteps_npt3}_{t_initial_npt3}-K_{t_final_npt3}-K_{tcoupl_npt3}_{ref_p_npt3}-bar_{pcoupl_npt3}" # Alterado de npt2 para npt3
        npt3_path = os.path.join(run_parent_dir, npt3_dir_name) # Alterado de npt2 para npt3
        os.makedirs(npt3_path, exist_ok=True)
        print(f"  Criado: '{npt3_path}'") # Alterado de npt2 para npt3

        # 4. NPT Stage 4 (Equilibration) - ANTERIORMENTE NPT Stage 3
        dt_npt4 = float(config_data['npt stage 4']['dt'])
        nsteps_npt4 = int(config_data['npt stage 4']['nsteps'])
        time_npt4 = (dt_npt4 * nsteps_npt4) / 1000 # Tempo em ns
        time_npt4_formatted = f"{time_npt4:.2f}".replace('.', 'p')
        t_initial_npt4 = config_data['npt stage 4']['t-initial']
        tcoupl_npt4 = config_data['npt stage 4']['tcoupl']
        ref_p_npt4 = config_data['npt stage 4']['ref-p']
        pcoupl_npt4 = config_data['npt stage 4']['pcoupl']
        npt4_dir_name = f"NPT_equilibration_{dt_npt4}_{time_npt4_formatted}ns_{nsteps_npt4}_{t_initial_npt4}-K_{tcoupl_npt4}_{ref_p_npt4}-bar_{pcoupl_npt4}"
        npt4_path = os.path.join(run_parent_dir, npt4_dir_name)
        os.makedirs(npt4_path, exist_ok=True)
        print(f"  Criado: '{npt4_path}'")

        # 5. NVT Stage 5 (Production) - ANTERIORMENTE NVT Stage 4
        dt_nvt5 = float(config_data['nvt stage 5']['dt'])
        nsteps_nvt5 = int(config_data['nvt stage 5']['nsteps'])
        time_nvt5 = (dt_nvt5 * nsteps_nvt5) / 1000 # Tempo em ns
        time_nvt5_formatted = f"{time_nvt5:.2f}".replace('.', 'p')
        t_initial_nvt5 = config_data['nvt stage 5']['t-initial']
        tcoupl_nvt5 = config_data['nvt stage 5']['tcoupl']
        nstxout_compressed_nvt5 = config_data['nvt stage 5']['nstxout-compressed']
        nvt5_dir_name = f"NVT_production_{dt_nvt5}_{time_nvt5_formatted}ns_{nsteps_nvt5}_{t_initial_nvt5}-K_{tcoupl_nvt5}_{nstxout_compressed_nvt5}-dump"
        nvt5_path = os.path.join(run_parent_dir, nvt5_dir_name)
        os.makedirs(nvt5_path, exist_ok=True)
        print(f"  Criado: '{nvt5_path}'")
        
        return min_dir_path # Retorna o caminho do diretório de minimização

    except KeyError as e:
        print(f"Erro: Variável '{e}' necessária para um dos diretórios de etapa não encontrada em user_variables.txt.")
        print("Por favor, verifique se todas as variáveis estão definidas e preenchidas nas seções corretas (ex: 'integrator' em [Energy Minimization]).")
        return None
    except ValueError as e:
        print(f"Erro de conversão de tipo: {e}. Verifique se os valores numéricos (dt, nsteps) no user_variables.txt são válidos.")
        return None
    except Exception as e:
        print(f"Um erro inesperado ocorreu ao criar os diretórios de etapa: {e}")
        return None


def _replace_mdp_parameters(template_content, params_dict):
    """
    Substitui parâmetros no conteúdo do template MDP com base no dicionário fornecido.
    Preserva linhas de comentário e outras linhas não modificadas.

    Args:
        template_content (str): O conteúdo completo do template.mdp.
        params_dict (dict): Dicionário de parâmetros a serem substituídos {param_name: new_value}.
                            As chaves do dicionário devem ser os nomes exatos dos parâmetros no MDP.

    Returns:
        str: O conteúdo do MDP com os parâmetros substituídos.
    """
    lines = template_content.splitlines()
    new_lines = []
    
    # Criar um conjunto de parâmetros que precisam ser substituídos para busca rápida
    params_to_replace = set(params_dict.keys())

    for line in lines:
        stripped_line = line.strip()
        # Ignorar linhas vazias ou comentários
        if not stripped_line or stripped_line.startswith(';') or stripped_line.startswith('#'):
            new_lines.append(line)
            continue

        replaced = False
        for param in params_to_replace:
            # Padrão regex para encontrar o parâmetro no início da linha,
            # seguido por um '=', e capturando o valor e qualquer comentário subsequente.
            # `re.escape(param)` é usado para tratar caracteres especiais no nome do parâmetro.
            # `\s*=\s*` permite espaços variados ao redor do '='.
            # `([^;]*)` captura o valor até um ';' ou o final da linha.
            # `(.*)` captura o restante da linha, incluindo comentários.
            pattern = rf"^{re.escape(param)}\s*=\s*([^;]*)(.*)$"
            match = re.match(pattern, line)
            
            if match:
                # O valor a ser inserido. Note que ConfigParser lê tudo como string.
                # Para 'annealing_time' e 'annealing_temp' que são listas de números,
                # o valor já virá como uma string com espaços, que é o esperado para o MDP.
                new_value = str(params_dict[param])
                comment_part = match.group(2) # O restante da linha (geralmente o comentário)

                # Tenta manter a formatação original da linha antes do '='
                eq_pos = line.find('=')
                if eq_pos != -1:
                    param_part = line[:eq_pos].strip()
                    leading_spaces_after_param = line[len(param_part):eq_pos]
                    new_line = f"{param_part}{leading_spaces_after_param}= {new_value}{comment_part}"
                else:
                    # Fallback caso não encontre '=', que não deveria acontecer para parâmetros MDP
                    new_line = f"{param:<25} = {new_value}{comment_part}" # Formatação padrão
                
                new_lines.append(new_line)
                replaced = True
                break # Move para a próxima linha do template, já que o parâmetro foi substituído
            
        if not replaced:
            new_lines.append(line) # Adicionar a linha original se nenhuma substituição ocorreu
            
    return "\n".join(new_lines)


def generate_mdp_files(main_sim_dir, config_data):
    """
    Gera os arquivos .mdp para cada etapa da simulação.
    """
    protocol_dir = os.path.join(main_sim_dir, "Protocol")
    
    # Localizar e ler o template.mdp
    # Assumimos que o MD_Simulation é um subdiretório do diretório onde o script é executado
    md_simulation_repo_path = os.path.join(os.getcwd(), "MD_Simulation")
    template_mdp_path = os.path.join(md_simulation_repo_path, "template.mdp")

    if not os.path.exists(template_mdp_path):
        print(f"Erro: O arquivo 'template.mdp' não foi encontrado em '{template_mdp_path}'.")
        print("Por favor, certifique-se de que o repositório 'MD_Simulation' foi clonado corretamente e 'template.mdp' está dentro dele.")
        return False

    try:
        with open(template_mdp_path, 'r') as f:
            template_content = f.read()
    except IOError as e:
        print(f"Erro ao ler o arquivo 'template.mdp': {e}")
        return False

    print("\nGerando arquivos .mdp para cada etapa da simulação...")

    # Dicionário de etapas e seus nomes de arquivo MDP, junto com as seções de user_variables
    # As chaves do dicionário de ConfigParser são convertidas para minúsculas.
    stages = {
        "Energy Minimization": {"file_name": "em.mdp", "section": "energy minimization"},
        "NVT 2 Stage":         {"file_name": "nvt_2.mdp", "section": "nvt 2 stage"},
        "NPT Stage 3":         {"file_name": "npt_3.mdp", "section": "npt stage 3"}, # Alterado de "NPT Stage 2" e "npt_2.mdp"
        "NPT Stage 4":         {"file_name": "npt_4.mdp", "section": "npt stage 4"},
        "NVT Stage 5":         {"file_name": "nvt_5.mdp", "section": "nvt stage 5"},
    }

    # Carrega os parâmetros gerais de MD (se a seção existir)
    md_general_params = {}
    if 'md_general' in config_data:
        md_general_params = {key: value for key, value in config_data['md_general'].items()}
    else:
        print("Aviso: Seção '[MD_General]' não encontrada em user_variables.txt. Parâmetros gerais de MD não serão aplicados automaticamente.")


    success = True
    for stage_display_name, info in stages.items():
        file_name = info["file_name"]
        section_name = info["section"]
        output_mdp_path = os.path.join(protocol_dir, file_name)

        if section_name not in config_data:
            print(f"Aviso: Seção `[{section_name.upper()}]` não encontrada em user_variables.txt para '{stage_display_name}'. Pulando a geração de '{file_name}'.")
            success = False
            continue

        # Começa com TODOS os parâmetros de MD_General como base
        stage_params_for_mdp = {key: value for key, value in md_general_params.items()}

        # Em seguida, sobrepõe os parâmetros específicos da etapa.
        # Isso garante que parâmetros específicos da etapa tenham precedência sobre os gerais,
        # mas que os gerais estejam presentes se não forem sobrescritos.
        specific_stage_config = config_data[section_name]
        for param, value in specific_stage_config.items():
            stage_params_for_mdp[param] = value # Isso adiciona ou sobrescreve

        # Gerar o conteúdo do MDP para esta etapa
        modified_content = _replace_mdp_parameters(template_content, stage_params_for_mdp)

        try:
            with open(output_mdp_path, 'w') as f:
                f.write(modified_content)
            print(f"  Gerado: '{output_mdp_path}' para '{stage_display_name}'.")
        except IOError as e:
            print(f"Erro ao escrever o arquivo '{output_mdp_path}': {e}")
            success = False
    
    return success

def create_subdirectories(parent_dir, config_data):
    """
    Cria os subdiretórios Force_Field, Protocol e RUN_{...} dentro do diretório pai.
    Retorna o caminho do diretório 'min_steep_10_500000' (ou None em caso de falha)
    e o caminho do diretório RUN principal.

    Args:
        parent_dir (str): O caminho para o diretório pai.
        config_data (configparser.ConfigParser): Objeto ConfigParser com as variáveis do usuário.
    Returns:
        tuple: (caminho_do_diretorio_min, caminho_do_diretorio_run_principal)
               Retorna (None, None) em caso de falha.
    """
    if not parent_dir:
        print("Não foi possível criar subdiretórios: o diretório pai não foi especificado.")
        return None, None # Indica falha

    print(f"\nCriando subdiretórios em '{parent_dir}'...")
    
    protocol_path = "" # Para armazenar o caminho da pasta Protocol
    run_parent_dir = None # Para armazenar o caminho da pasta RUN_...
    min_dir_path_return = None # Para armazenar o caminho da pasta de minimização

    fixed_subdirs = ["Force_Field", "Protocol"]

    for subdir in fixed_subdirs:
        path = os.path.join(parent_dir, subdir)
        try:
            os.makedirs(path, exist_ok=True)
            print(f"  Criado: '{path}'")
            if subdir == "Protocol":
                protocol_path = path # Guarda o caminho da pasta Protocol
        except OSError as e:
            print(f"Erro ao criar o diretório '{path}': {e}")
            return None, None

    try:
        system = config_data['general']['system']
        salt_concentration = config_data['general'].get('salt_concentration', 'NoSalt') 
        t_target = config_data['general'].get('t_target', 'N_A') 

        run_dir_name = f"RUN_{system}_{salt_concentration}_{t_target}"
        run_path = os.path.join(parent_dir, run_dir_name)
        os.makedirs(run_path, exist_ok=True)
        print(f"  Criado: '{run_path}'")
        run_parent_dir = run_path

    except KeyError as e:
        print(f"Erro: Variável '{e}' necessária para o nome do diretório RUN não encontrada ou vazia.")
        print("Por favor, verifique se 'system', 'Salt_concentration' e 'T_target' estão definidos na seção [General] de user_variables.txt.")
        return None, None
    except Exception as e:
        print(f"Um erro inesperado ocorreu ao criar o diretório RUN: {e}")
        return None, None
    
    # Chama a função para criar os subdiretórios de etapa dentro de RUN
    if run_parent_dir:
        min_dir_path_return = create_run_subdirectories_actual(run_parent_dir, config_data)
        if not min_dir_path_return: # Se create_run_subdirectories_actual retornou None
            print("Não foi possível criar todos os subdiretórios de etapa de RUN. Encerrando.")
            return None, None # Indica falha
    
    return min_dir_path_return, run_parent_dir # Retorna o caminho do diretório de minimização e do RUN principal

def create_notes_references_directory(main_sim_dir):
    """
    Cria o diretório 'Notes_References' dentro do diretório principal da simulação.

    Args:
        main_sim_dir (str): O caminho para o diretório principal da simulação.

    Returns:
        bool: True se o diretório foi criado ou já existia, False em caso de erro.
    """
    if not main_sim_dir:
        print("Não foi possível criar o diretório 'Notes_References': o diretório principal da simulação não foi especificado.")
        return False

    notes_dir_path = os.path.join(main_sim_dir, "Notes_References")

    try:
        os.makedirs(notes_dir_path, exist_ok=True)
        print(f"  Criado: '{notes_dir_path}'")
        return True
    except OSError as e:
        print(f"Erro ao criar o diretório 'Notes_References' em '{notes_dir_path}': {e}")
        return False
    except Exception as e:
        print(f"Um erro inesperado ocorreu ao criar o diretório 'Notes_References': {e}")
        return False

def copy_cpu_template_job(run_target_dir):
    """
    Copia o arquivo 'CPU-template.job' do repositório clonado para o
    diretório de destino (o diretório RUN_*), renomeando-o para 'cpu.job'.

    Args:
        run_target_dir (str): O caminho para o diretório RUN_*.
    Returns:
        bool: True se a cópia foi bem-sucedida, False caso contrário.
    """
    if not run_target_dir:
        print("Não foi possível copiar 'CPU-template.job': o diretório RUN_* não foi especificado.")
        return False

    md_simulation_repo_path = os.path.join(os.getcwd(), "MD_Simulation")
    source_cpu_template_path = os.path.join(md_simulation_repo_path, "CPU-template.job")
    output_cpu_job_path = os.path.join(run_target_dir, "cpu.job")

    if not os.path.exists(source_cpu_template_path):
        print(f"Erro: O arquivo 'CPU-template.job' não foi encontrado em '{source_cpu_template_path}'.")
        print("Certifique-se de que o repositório 'MD_Simulation' foi clonado corretamente e o arquivo existe.")
        return False

    try:
        shutil.copy2(source_cpu_template_path, output_cpu_job_path)
        print(f"\nCopiado 'CPU-template.job' para '{output_cpu_job_path}' como 'cpu.job'.")
        return True
    except Exception as e:
        print(f"Erro inesperado ao copiar o arquivo 'CPU-template.job': {e}")
        return False

def copy_pdb_files(config_data, destination_dir):
    """
    Copia arquivos .pdb das pastas de origem para o diretório de destino.
    Os nomes dos arquivos .pdb são baseados em 'molecule1', 'molecule2', 'molecule3'.

    Args:
        config_data (configparser.ConfigParser): Objeto ConfigParser com as variáveis do usuário.
        destination_dir (str): O caminho para o diretório de destino (e.g., min_steep_10_500000).
    Returns:
        bool: True se todos os arquivos foram processados (copiados ou avisados), False se houver um erro grave.
    """
    if 'paths' not in config_data:
        print("Aviso: Seção '[Paths]' não encontrada em user_variables.txt. Pulando a cópia de arquivos PDB.")
        return False

    pdb_source_folder = config_data['paths'].get('original_pdb_folder', '').strip()
    if not pdb_source_folder:
        print("Aviso: 'Original_PDB_Folder' não está definido ou está vazio na seção [Paths]. Pulando a cópia de arquivos PDB.")
        return False
    
    if not os.path.isdir(pdb_source_folder):
        print(f"Erro: O diretório de origem PDB '{pdb_source_folder}' não existe ou não é um diretório válido. Verifique o caminho em user_variables.txt.")
        return False

    molecules = []
    # Usar get com valor padrão para evitar KeyError se a variável não estiver presente
    for i in range(1, 4): # Checa molecule1, molecule2, molecule3
        mol_name = config_data['general'].get(f'molecule{i}', '').strip()
        if mol_name:
            molecules.append(mol_name)
    
    if not molecules:
        print("Aviso: Nenhuma molécula (molecule1, molecule2, molecule3) definida na seção [General]. Nenhuns arquivos PDB para copiar.")
        return False

    print(f"\nCopiando arquivos .pdb de '{pdb_source_folder}' para '{destination_dir}'...")
    all_copied = True
    for mol in molecules:
        source_pdb_path = os.path.join(pdb_source_folder, f"{mol}.pdb")
        dest_pdb_path = os.path.join(destination_dir, f"{mol}.pdb")
        if os.path.exists(source_pdb_path):
            try:
                shutil.copy2(source_pdb_path, dest_pdb_path) # copy2 preserva metadados
                print(f"  Copiado: '{source_pdb_path}' para '{dest_pdb_path}'")
            except Exception as e:
                print(f"Erro ao copiar '{source_pdb_path}' para '{dest_pdb_path}': {e}")
                all_copied = False
        else:
            print(f"Aviso: Arquivo PDB '{source_pdb_path}' para a molécula '{mol}' não encontrado. Pulando a cópia.")
            all_copied = False
    
    return all_copied

def copy_packmol_template_only(destination_dir):
    """
    Copia o arquivo 'packmol-template.inp' do repositório clonado para o
    diretório de destino (o diretório de minimização), renomeando-o para 'packmol.inp'.
    Esta função NÃO realiza substituições ou inserções de blocos de texto no arquivo.

    Args:
        destination_dir (str): O caminho para o diretório de destino (e.g., min_steep_10_500000).
    Returns:
        bool: True se a cópia foi bem-sucedida, False caso contrário.
    """
    md_simulation_repo_path = os.path.join(os.getcwd(), "MD_Simulation")
    template_packmol_path = os.path.join(md_simulation_repo_path, "packmol-template.inp")
    output_packmol_path = os.path.join(destination_dir, "packmol.inp")

    if not os.path.exists(template_packmol_path):
        print(f"Erro: O arquivo 'packmol-template.inp' não foi encontrado em '{template_packmol_path}'.")
        print("Certifique-se de que o repositório 'MD_Simulation' foi clonado corretamente e o arquivo existe.")
        return False

    try:
        # Apenas copia o arquivo, sem nenhuma alteração no conteúdo
        shutil.copy2(template_packmol_path, output_packmol_path)
        print(f"\nCopiado 'packmol-template.inp' para '{output_packmol_path}' sem modificações adicionais.")
        return True

    except Exception as e:
        print(f"Erro inesperado ao copiar o arquivo Packmol: {e}")
        return False

def process_topology_template(config_data, destination_dir):
    """
    Copia o arquivo 'topology-template.top' do repositório clonado para o
    diretório de destino (o diretório de minimização), renomeando-o para 'topology.top',
    e substitui variáveis dentro de chaves {} pelos seus valores do user_variables.txt.

    Args:
        config_data (configparser.ConfigParser): Objeto ConfigParser com as variáveis do usuário.
        destination_dir (str): O caminho para o diretório de destino (e.g., min_steep_10_500000).
    Returns:
        bool: True se o processamento foi bem-sucedido, False caso contrário.
    """
    md_simulation_repo_path = os.path.join(os.getcwd(), "MD_Simulation")
    template_top_path = os.path.join(md_simulation_repo_path, "topology-template.top")
    output_top_path = os.path.join(destination_dir, "topology.top")

    if not os.path.exists(template_top_path):
        print(f"Erro: O arquivo 'topology-template.top' não foi encontrado em '{template_top_path}'.")
        print("Certifique-se de que o repositório 'MD_Simulation' foi clonado corretamente e o arquivo existe.")
        return False

    try:
        # 1. Ler o conteúdo do template
        with open(template_top_path, 'r') as f:
            template_content = f.read()

        # Criar um dicionário de todos os parâmetros relevantes do config_data para fácil lookup
        # ConfigParser já converte chaves para minúsculas
        all_params = {}
        for section in config_data.sections():
            for key, value in config_data[section].items():
                all_params[key] = value.strip() # Remove espaços em branco do valor

        def replace_match(match):
            var_name = match.group(1) # O nome da variável dentro das chaves (e.g., 'system', 'molecule1')
            
            # Tenta obter o valor do dicionário de todos os parâmetros
            value = all_params.get(var_name.lower(), None) # .lower() para corresponder às chaves do ConfigParser
            
            if value is not None:
                if value == "": # Se a variável está vazia no user_variables.txt
                    print(f"Aviso: Variável TOPOLOGY '{var_name}' está vazia em user_variables.txt. Substituindo por string vazia.")
                    return ""
                return value
            else:
                print(f"Aviso: Variável TOPOLOGY '{var_name}' não encontrada em user_variables.txt. Deixando o placeholder original.")
                return match.group(0) # Retorna a string original (ex: "{system}") se não encontrar

        # Usa re.sub com uma função de substituição para lidar com cada ocorrência de {VARIAVEL}
        modified_content = re.sub(r"\{([a-zA-Z0-9_]+)\}", replace_match, template_content)

        # 3. Escrever o conteúdo modificado para o arquivo de saída
        with open(output_top_path, 'w') as f:
            f.write(modified_content)
            
        print(f"\nProcessado e copiado 'topology-template.top' para '{output_top_path}' com variáveis substituídas.")
        return True

    except Exception as e:
        print(f"Erro inesperado ao processar o arquivo de topologia: {e}")
        return False

def copy_itp_files(config_data, force_field_dir):
    """
    Copia arquivos .itp da pasta de origem (Original_FF_Folder) para o diretório
    Force_Field, com base nos nomes das moléculas definidos em user_variables.txt
    e em nomes de arquivos ITP comuns.

    Args:
        config_data (configparser.ConfigParser): Objeto ConfigParser com as variáveis do usuário.
        force_field_dir (str): O caminho para o diretório Force_Field criado pelo script.
    Returns:
        bool: True se todos os arquivos relevantes foram copiados ou avisados, False se houver um erro grave.
    """
    if 'paths' not in config_data:
        print("Aviso: Seção '[Paths]' não encontrada em user_variables.txt. Pulando a cópia de arquivos ITP.")
        return False

    ff_source_folder = config_data['paths'].get('original_ff_folder', '').strip()
    if not ff_source_folder:
        print("Aviso: 'Original_FF_Folder' não está definido ou está vazio na seção [Paths]. Pulando a cópia de arquivos ITP.")
        return False
    
    if not os.path.isdir(ff_source_folder):
        print(f"Erro: O diretório de origem do Force Field '{ff_source_folder}' não existe ou não é um diretório válido. Verifique o caminho em user_variables.txt.")
        return False

    # Coletar todos os nomes de moléculas definidos na seção [general]
    molecule_names = []
    for key, value in config_data['general'].items():
        if key.startswith('molecule') and value.strip(): # Garante que a molécula não é vazia
            molecule_names.append(value.strip())

    # Lista de arquivos ITP comuns a serem procurados e copiados
    # Inclui o nome dos sistemas de force field (ex: amber99sb-ildn.ff)
    common_itp_files = [
        "ions.itp",
        "tip3p.itp",
        config_data['general'].get('force_field', '').strip() + ".ff/forcefield.itp", # Ex: amber99sb-ildn.ff/forcefield.itp
        config_data['general'].get('force_field', '').strip() + ".ff/ions.itp", # Para caso ions.itp esteja dentro da pasta ff
        config_data['general'].get('force_field', '').strip() + ".ff/tip3p.itp" # Para caso tip3p.itp esteja dentro da pasta ff
    ]
    # Adicionar itp para cada molécula
    for mol in molecule_names:
        common_itp_files.append(f"{mol}.itp")
    
    # Adicionar o arquivo de topologia principal do force field (ex: amber99sb-ildn.ff/forcefield.itp)
    # Certificar-se que o nome do FF está presente em user_variables.txt
    force_field_name = config_data['general'].get('force_field', '').strip()
    if force_field_name:
        common_itp_files.append(os.path.join(f"{force_field_name}.ff", "forcefield.itp"))
        # Adicionar também o arquivo de empacotamento específico do force field se houver (ex: water.itp)
        # Este é um exemplo, pode precisar de mais lógica se houver nomes variáveis
        # Por exemplo, se a água não for tip3p, mas algum outro modelo.
        # Por agora, assumimos que 'tip3p.itp' cobre a água.

    print(f"\nCopiando arquivos .itp de '{ff_source_folder}' para '{force_field_dir}'...")
    all_copied = True
    for itp_file in common_itp_files:
        source_itp_path = os.path.join(ff_source_folder, itp_file)
        dest_itp_path = os.path.join(force_field_dir, os.path.basename(itp_file)) # Copia para o diretório FF
        
        # Para lidar com os arquivos dentro de subpastas (como forcefield.itp dentro de amber99sb-ildn.ff)
        # Precisamos garantir que a estrutura de subpastas dentro de Force_Field seja criada se necessário.
        # No entanto, a instrução é copiar para o 'Force_Field' principal.
        # Vamos ajustar para copiar direto para 'Force_Field' e assumir que o 'include' no .top resolve.
        # Se for o caso de amber99sb-ildn.ff/forcefield.itp, só copia o forcefield.itp para Force_Field/forcefield.itp
        
        if os.path.exists(source_itp_path):
            try:
                # Criar o diretório de destino se não existir (se for um subdiretório dentro de force_field_dir)
                # No momento, estamos copiando tudo diretamente para force_field_dir, então não precisa disso.
                # Se fosse para manter a estrutura de subdiretórios, precisaríamos:
                # os.makedirs(os.path.dirname(dest_itp_path), exist_ok=True)
                shutil.copy2(source_itp_path, dest_itp_path)
                print(f"  Copiado: '{source_itp_path}' para '{dest_itp_path}'")
            except Exception as e:
                print(f"Erro ao copiar '{source_itp_path}' para '{dest_itp_path}': {e}")
                all_copied = False
        else:
            # Não é um erro se um ITP comum não for encontrado (exceto para as moléculas)
            # Apenas avisa se for um arquivo ITP importante (molécula ou forcefield.itp)
            if any(m_name in itp_file for m_name in molecule_names) or "forcefield.itp" in itp_file:
                 print(f"Aviso: Arquivo ITP crítico '{source_itp_path}' não encontrado. Pode ser um erro na topologia.")
            else:
                 print(f"Aviso: Arquivo ITP '{source_itp_path}' não encontrado. Pulando a cópia.")
            all_copied = False
    
    return all_copied
    
# --- Função Principal ---
def main():
    repo_url = "git@github.com:tuananlourenco/MD_Simulation.git"
    user_variables_file = "user_variables.txt"

    # Clonar o repositório
    if not clone_repository(repo_url):
        print("Falha na clonagem do repositório ou o diretório já existe. Verifique as mensagens acima. Encerrando.")
        return

    # Ler user_variables.txt
    config_data = read_user_variables(user_variables_file)
    if config_data is None:
        print("Falha ao ler o arquivo de variáveis do usuário. Encerrando.")
        return

    # Criar o diretório principal de simulação (Project_...)
    main_simulation_directory = create_main_simulation_directory(config_data)
    if not main_simulation_directory:
        print("Falha ao criar o diretório principal de simulação. Encerrando.")
        return

    # Criar subdiretórios (Force_Field, Protocol, RUN_...)
    # min_dir_path: Caminho para o diretório de minimização (ex: min_steep_10_500000)
    # run_parent_dir: Caminho para o diretório RUN_...
    min_dir_path, run_parent_dir = create_subdirectories(main_simulation_directory, config_data)
    if min_dir_path is None or run_parent_dir is None:
        print("Falha ao criar subdiretórios essenciais. Encerrando.")
        return
    
    # Criar diretório Notes_References
    if not create_notes_references_directory(main_simulation_directory):
        print("Falha ao criar o diretório Notes_References. Continuando, mas pode causar problemas.")
        # Não encerra o script, pois este não é um diretório crítico para a simulação em si

    # Gerar arquivos .mdp
    if not generate_mdp_files(main_simulation_directory, config_data):
        print("Alguns arquivos MDP não puderam ser gerados. Verifique os avisos acima. Continuando, mas pode causar problemas.")
        # Não encerra, mas avisa

    # Copiar CPU-template.job para cpu.job no diretório RUN_...
    if not copy_cpu_template_job(run_parent_dir):
        print("Falha ao copiar 'CPU-template.job'. Continuando, mas o script de job pode estar faltando.")
        # Não encerra, mas avisa

    # Copiar arquivos PDB para o diretório de minimização
    if not copy_pdb_files(config_data, min_dir_path):
        print("Alguns arquivos PDB não puderam ser copiados. Verifique os avisos acima. Continuando, mas isso pode ser crítico para a simulação.")
        # Não encerra, mas avisa. Isso PODE ser crítico, dependendo se as moléculas faltantes são essenciais.

    # Copiar packmol-template.inp para packmol.inp no diretório de minimização
    if not copy_packmol_template_only(min_dir_path):
        print("Falha ao copiar o arquivo Packmol template. Continuando, mas o Packmol pode não funcionar.")
        # Não encerra, mas avisa

    # Processar e copiar topology-template.top para topology.top no diretório de minimização
    if not process_topology_template(config_data, min_dir_path):
        print("Falha ao processar e copiar o arquivo de topologia. Continuando, mas isso é crítico para a simulação.")
        # Não encerra, mas avisa. Isso é CRÍTICO.

    # Copiar arquivos ITP para o diretório Force_Field
    force_field_dir = os.path.join(main_simulation_directory, "Force_Field")
    if not copy_itp_files(config_data, force_field_dir):
        print("Alguns arquivos ITP não puderam ser copiados. Verifique os avisos acima. Isso é crítico para a simulação.")
        # Não encerra, mas avisa. Isso é CRÍTICO.
    
    print("\nConfiguração inicial da simulação concluída. Verifique as mensagens acima para quaisquer avisos ou erros.")

if __name__ == "__main__":
    main()
