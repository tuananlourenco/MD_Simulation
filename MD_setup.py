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
    target_path = os.path.join(os.getcwd(), clone_dir)

    if os.path.exists(target_path) and os.listdir(target_path):
        print(f"O repositório '{repo_url}' parece já ter sido clonado em '{target_path}'. Pulando a clonagem.")
        return True

    print(f"Clonando o repositório '{repo_url}' para '{target_path}'...")
    try:
        # Adicione o argumento '--depth 1' para clonar apenas o histórico mais recente, se o repositório for grande.
        # Caso o acesso seja via SSH (git@github.com), a autenticação deve estar configurada (chaves SSH).
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

        # 3. NPT Stage 2 (Annealing)
        dt_npt2 = float(config_data['npt stage 2']['dt'])
        nsteps_npt2 = int(config_data['npt stage 2']['nsteps'])
        time_npt2 = (dt_npt2 * nsteps_npt2) / 1000 # Tempo em ns
        time_npt2_formatted = f"{time_npt2:.2f}".replace('.', 'p')
        t_initial_npt2 = config_data['npt stage 2']['t-initial']
        t_final_npt2 = config_data['npt stage 2']['t-final']
        tcoupl_npt2 = config_data['npt stage 2']['tcoupl']
        ref_p_npt2 = config_data['npt stage 2']['ref-p']
        pcoupl_npt2 = config_data['npt stage 2']['pcoupl']
        npt2_dir_name = f"NPT_annealing_{dt_npt2}_{time_npt2_formatted}ns_{nsteps_npt2}_{t_initial_npt2}-K_{t_final_npt2}-K_{tcoupl_npt2}_{ref_p_npt2}-bar_{pcoupl_npt2}"
        npt2_path = os.path.join(run_parent_dir, npt2_dir_name)
        os.makedirs(npt2_path, exist_ok=True)
        print(f"  Criado: '{npt2_path}'")

        # 4. NPT Stage 3 (Equilibration)
        dt_npt3 = float(config_data['npt stage 3']['dt'])
        nsteps_npt3 = int(config_data['npt stage 3']['nsteps'])
        time_npt3 = (dt_npt3 * nsteps_npt3) / 1000 # Tempo em ns
        time_npt3_formatted = f"{time_npt3:.2f}".replace('.', 'p')
        t_initial_npt3 = config_data['npt stage 3']['t-initial']
        tcoupl_npt3 = config_data['npt stage 3']['tcoupl']
        ref_p_npt3 = config_data['npt stage 3']['ref-p']
        pcoupl_npt3 = config_data['npt stage 3']['pcoupl']
        npt3_dir_name = f"NPT_equilibration_{dt_npt3}_{time_npt3_formatted}ns_{nsteps_npt3}_{t_initial_npt3}-K_{tcoupl_npt3}_{ref_p_npt3}-bar_{pcoupl_npt3}"
        npt3_path = os.path.join(run_parent_dir, npt3_dir_name)
        os.makedirs(npt3_path, exist_ok=True)
        print(f"  Criado: '{npt3_path}'")

        # 5. NVT Stage 4 (Production)
        dt_nvt4 = float(config_data['nvt stage 4']['dt'])
        nsteps_nvt4 = int(config_data['nvt stage 4']['nsteps'])
        time_nvt4 = (dt_nvt4 * nsteps_nvt4) / 1000 # Tempo em ns
        time_nvt4_formatted = f"{time_nvt4:.2f}".replace('.', 'p')
        t_initial_nvt4 = config_data['nvt stage 4']['t-initial']
        tcoupl_nvt4 = config_data['nvt stage 4']['tcoupl']
        nstxout_compressed_nvt4 = config_data['nvt stage 4']['nstxout-compressed']
        nvt4_dir_name = f"NVT_production_{dt_nvt4}_{time_nvt4_formatted}ns_{nsteps_nvt4}_{t_initial_nvt4}-K_{tcoupl_nvt4}_{nstxout_compressed_nvt4}-dump"
        nvt4_path = os.path.join(run_parent_dir, nvt4_dir_name)
        os.makedirs(nvt4_path, exist_ok=True)
        print(f"  Criado: '{nvt4_path}'")
        
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
        "NPT Stage 2":         {"file_name": "npt_2.mdp", "section": "npt stage 2"},
        "NPT Stage 3":         {"file_name": "npt_3.mdp", "section": "npt stage 3"},
        "NVT Stage 4":         {"file_name": "nvt_4.mdp", "section": "nvt stage 4"},
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
        if key.startswith('molecule') and value.strip():
            molecule_names.append(value.strip())
    
    if not molecule_names:
        print("Aviso: Nenhuma molécula (e.g., molecule1, molecule2) definida e preenchida na seção [General]. Arquivos ITP baseados em nomes de moléculas não serão copiados.")
        # O script continua para verificar ITPs comuns mesmo sem moléculas definidas
    
    print(f"\nCopiando arquivos .itp de '{ff_source_folder}' para '{force_field_dir}'...")
    all_copied = True
    
    # Lista de arquivos ITP no diretório de origem
    itp_files_in_source = [f for f in os.listdir(ff_source_folder) if f.lower().endswith('.itp')]

    # Lista de nomes comuns de arquivos ITP que devem ser copiados independentemente do nome da molécula
    common_itps = ["forcefield.itp", "ions.itp", "water.itp"] # Adicione mais conforme necessário

    for filename in itp_files_in_source:
        source_itp_path = os.path.join(ff_source_folder, filename)
        dest_itp_path = os.path.join(force_field_dir, filename)
        
        should_copy = False
        
        # 1. Verifica se é um arquivo ITP comum
        if filename.lower() in common_itps:
            should_copy = True
        else:
            # 2. Verifica se o nome do arquivo ITP contém o nome de alguma das moléculas definidas
            for mol_name in molecule_names:
                if mol_name.lower() in filename.lower(): # Comparação case-insensitive
                    should_copy = True
                    break # Encontrou uma correspondência para este arquivo, não precisa verificar outras moléculas
        
        if should_copy:
            if os.path.exists(source_itp_path):
                try:
                    shutil.copy2(source_itp_path, dest_itp_path) # copy2 preserva metadados
                    print(f"  Copiado: '{source_itp_path}' para '{dest_itp_path}'")
                except Exception as e:
                    print(f"Erro ao copiar '{source_itp_path}' para '{dest_itp_path}': {e}")
                    all_copied = False
            else:
                print(f"Aviso: Arquivo ITP '{source_itp_path}' não encontrado. Pulando a cópia.")
                all_copied = False
        else:
            # Mensagem informativa para arquivos ITP que não são relevantes para cópia
            print(f"Informação: Arquivo ITP '{filename}' não corresponde a nenhuma molécula definida ou arquivo ITP comum. Pulando.")

    return all_copied


if __name__ == "__main__":
    # URL do repositório GitHub para clonar
    github_repo_url = "git@github.com:tuananlourenco/MD_Simulation.git" 
    
    # Nome do arquivo de variáveis do usuário
    user_variables_file = "user_variables.txt"

    # 1. Clonar o repositório
    if not clone_repository(github_repo_url):
        print("\nFalha na clonagem do repositório ou o diretório já existe. Verifique as mensagens acima. Encerrando.")
        exit(1)

    # 2. Ler as variáveis do usuário
    config_data = read_user_variables(user_variables_file)
    if not config_data:
        print("\nNão foi possível carregar as variáveis do usuário. Encerrando.")
        exit(1)

    # 3. Criar o diretório principal da simulação
    main_sim_dir = create_main_simulation_directory(config_data)
    if not main_sim_dir:
        print("\nNão foi possível criar o diretório principal de simulação. Encerrando.")
        exit(1)

    # 4. Criar os subdiretórios (Force_Field, Protocol, RUN_...)
    # create_subdirectories agora retorna o caminho do diretório 'min_steep_10_500000'
    min_em_dir_path, main_run_dir_path = create_subdirectories(main_sim_dir, config_data)
    if not min_em_dir_path: # Se min_em_dir_path for None, algo falhou na criação dos diretórios
        print("\nNão foi possível criar todos os subdiretórios necessários. Encerrando.")
        exit(1)

    # Definir o caminho para o diretório Force_Field (criado em create_subdirectories)
    force_field_dir_path = os.path.join(main_sim_dir, "Force_Field")

    # 5. Gerar os arquivos .mdp
    if not generate_mdp_files(main_sim_dir, config_data):
        print("\nHouve erros na geração de alguns arquivos .mdp. Verifique as mensagens acima.")
        # O script continua para outras etapas mesmo com erros no MDP, mas você pode mudar para exit(1) se quiser parar.

    # 6. Copiar arquivos PDB para o diretório de minimização
    if not copy_pdb_files(config_data, min_em_dir_path):
        print("\nHouve avisos ou erros durante a cópia dos arquivos PDB. Verifique as mensagens acima.")
        # O script continua mesmo com erros na cópia de PDBs.

    # 7. Copiar o template do Packmol (apenas cópia, sem modificações de conteúdo)
    if not copy_packmol_template_only(min_em_dir_path):
        print("\nHouve erros ao copiar o template do Packmol. Verifique as mensagens acima.")
        # O script continua mesmo com erros na cópia do Packmol.

    # 8. Processar o template da Topologia (substitui variáveis)
    if not process_topology_template(config_data, min_em_dir_path):
        print("\nHouve erros ao processar o template da Topologia. Verifique as mensagens acima.")
        # O script continua mesmo com erros na topologia.

    # 9. Copiar arquivos ITP para o diretório Force_Field
    if not copy_itp_files(config_data, force_field_dir_path):
        print("\nHouve avisos ou erros ao copiar arquivos ITP. Verifique as mensagens acima.")
        # O script continua mesmo com erros na cópia de ITPs.


    print(f"\nProcesso completo! Estrutura de diretórios, arquivos .mdp, PDBs, Packmol.inp, topology.top e ITPs (se encontrados) criados/atualizados com sucesso em: {main_sim_dir}")

