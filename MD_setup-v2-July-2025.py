import sys
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

class UnitStrippingConfigParser(configparser.ConfigParser):
    """
    Uma classe personalizada que herda de ConfigParser.
    Ela modifica o método 'get' para remover automaticamente unidades
    e comentários entre chaves (ex: {K}, {nm}) dos valores lidos.
    """
    def get(self, section, option, *, raw=False, vars=None, fallback=configparser._UNSET):
        # Primeiro, obtemos o valor original usando o método da classe pai
        value = super().get(section, option, raw=raw, vars=vars, fallback=fallback)

        # Se o valor for uma string (o caso mais comum), removemos a parte com a unidade
        if isinstance(value, str):
            # Usamos split('{')[0] para pegar tudo antes do '{' e .strip() para remover espaços
            return value.split('{')[0].strip()

        # Se não for uma string (ou se for o valor de fallback), retornamos como está
        return value

def read_user_variables(file_path):
    """
    Lê o arquivo user_variables.txt e armazena as variáveis em um objeto 
    UnitStrippingConfigParser, que remove automaticamente as unidades em {}.
    Ignora comentários (linhas começando com # ou ;) e linhas vazias.

    Args:
        file_path (str): O caminho para o arquivo user_variables.txt.

    Returns:
        configparser.ConfigParser: Um objeto UnitStrippingConfigParser contendo as variáveis limpas.
                                   Retorna None se o arquivo não for encontrado ou houver erro.
    """
    # --- ALTERAÇÃO PRINCIPAL AQUI ---
    # Usamos nossa nova classe em vez do configparser.ConfigParser padrão
    user_vars = UnitStrippingConfigParser(allow_no_value=True)
    
    try:
        read_files = user_vars.read(file_path)
        if not read_files:
            print(f"Aviso: O arquivo '{file_path}' foi encontrado, mas nenhuma seção válida foi lida.")
            return None

        print(f"Variáveis lidas de '{file_path}' com sucesso! As unidades em {{}} serão ignoradas.")
        print(f"Seções encontradas pelo ConfigParser: {user_vars.sections()}")

        if 'general' not in user_vars.sections():
            print("\nERRO DE DEPURACÃO: A seção '[general]' não foi encontrada pelo ConfigParser.")
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
    Cria os subdiretórios específicos de etapas da simulação dentro do diretório RUN
    com a nova nomenclatura numerada e detalhada.
    """
    if not run_parent_dir:
        print("Não foi possível criar diretórios de etapas: o diretório RUN não foi especificado.")
        return None

    print(f"\nCriando diretórios de etapas da simulação em '{run_parent_dir}'...")

    min_dir_path = None

    try:
        # --- 1. Minimização ---
        section = '1-min'
        integrator_em = config_data.get(section, 'integrator', fallback='NA')
        emtol_em = config_data.get(section, 'emtol', fallback='NA')
        emtol_em_str = emtol_em.replace('.', 'p')
        nsteps_em = config_data.get(section, 'nsteps', fallback='NA')
        emstep_em = config_data.get(section, 'emstep', fallback='NA')
        emstep_em_str = emstep_em.replace('.', 'p') # CORREÇÃO AQUI
        em_dir_name = f"1-min_integrator-{integrator_em}_emtol-{emtol_em_str}_nsteps-{nsteps_em}_emstep-{emstep_em_str}"
        em_path = os.path.join(run_parent_dir, em_dir_name)
        os.makedirs(em_path, exist_ok=True)
        print(f"  Criado: '{em_path}'")
        min_dir_path = em_path
        # --- 2. NVT Thermalization ---
        section = '2-NVT_thermalization'
        integrator_nvt2 = config_data.get(section, 'integrator', fallback='NA')
        dt_nvt2 = float(config_data.get(section, 'dt', fallback=0))
        dt_nvt2_str = str(dt_nvt2).replace('.', 'p')
        nsteps_nvt2 = int(config_data.get(section, 'nsteps', fallback=0))
        time_nvt2_ps = int(dt_nvt2 * nsteps_nvt2)
        ref_t_nvt2 = config_data.get(section, 'ref-t', fallback='NA')
        tcoupl_nvt2 = config_data.get(section, 'tcoupl', fallback='NA')
        nvt2_dir_name = f"2-NVT_thermalization_integrator-{integrator_nvt2}_dt_{dt_nvt2_str}_time_{time_nvt2_ps}ps_T-{ref_t_nvt2}_{tcoupl_nvt2}"
        nvt2_path = os.path.join(run_parent_dir, nvt2_dir_name)
        os.makedirs(nvt2_path, exist_ok=True)
        print(f"  Criado: '{nvt2_path}'")

        # --- 3. NPT Ergodicity ---
        section = '3-NPT_ergodicity'
        integrator_npt3 = config_data.get(section, 'integrator', fallback='NA')
        dt_npt3 = float(config_data.get(section, 'dt', fallback=0))
        dt_npt3_str = str(dt_npt3).replace('.', 'p')
        nsteps_npt3 = int(config_data.get(section, 'nsteps', fallback=0))
        time_npt3_ps = int(dt_npt3 * nsteps_npt3)
        t_initial_npt3 = config_data.get(section, 't-initial', fallback='NA')
        t_final_npt3 = config_data.get(section, 't-final', fallback='NA')
        ref_p_npt3 = config_data.get(section, 'ref-p', fallback='NA')
        ref_p_npt3_str = ref_p_npt3.replace('.', 'p')
        tcoupl_npt3 = config_data.get(section, 'tcoupl', fallback='NA')
        pcoupl_npt3 = config_data.get(section, 'pcoupl', fallback='NA')
        npt3_dir_name = f"3-NPT_ergodicity_integrator-{integrator_npt3}_dt_{dt_npt3_str}_time_{time_npt3_ps}ps_T-{t_initial_npt3}-{t_final_npt3}-K_P-{ref_p_npt3_str}-bar_{tcoupl_npt3}_{pcoupl_npt3}"
        npt3_path = os.path.join(run_parent_dir, npt3_dir_name)
        os.makedirs(npt3_path, exist_ok=True)
        print(f"  Criado: '{npt3_path}'")

        # --- 4. NPT Equilibration ---
        section = '4-NPT_equilibration'
        integrator_npt4 = config_data.get(section, 'integrator', fallback='NA')
        dt_npt4 = float(config_data.get(section, 'dt', fallback=0))
        dt_npt4_str = str(dt_npt4).replace('.', 'p')
        nsteps_npt4 = int(config_data.get(section, 'nsteps', fallback=0))
        time_npt4_ps = int(dt_npt4 * nsteps_npt4)
        ref_t_npt4 = config_data.get(section, 'ref-t', fallback='NA')
        ref_p_npt4 = config_data.get(section, 'ref-p', fallback='NA')
        ref_p_npt4_str = ref_p_npt4.replace('.', 'p')
        tcoupl_npt4 = config_data.get(section, 'tcoupl', fallback='NA')
        pcoupl_npt4 = config_data.get(section, 'pcoupl', fallback='NA')
        npt4_dir_name = f"4-NPT_equilibration_integrator-{integrator_npt4}_dt_{dt_npt4_str}_time_{time_npt4_ps}ps_T-{ref_t_npt4}-K_P-{ref_p_npt4_str}-bar_{tcoupl_npt4}_{pcoupl_npt4}"
        npt4_path = os.path.join(run_parent_dir, npt4_dir_name)
        os.makedirs(npt4_path, exist_ok=True)
        print(f"  Criado: '{npt4_path}'")
        
        # --- NOVA ETAPA: 6. NVT Re-equilibrium ---
        section = '6-NVT_re-equilibrium'
        if section in config_data:
            integrator_nvt6 = config_data.get(section, 'integrator', fallback='NA')
            dt_nvt6 = float(config_data.get(section, 'dt', fallback=0))
            dt_nvt6_str = str(dt_nvt6).replace('.', 'p')
            nsteps_nvt6 = int(config_data.get(section, 'nsteps', fallback=0))
            time_nvt6_ps = int(dt_nvt6 * nsteps_nvt6)
            ref_t_nvt6 = config_data.get(section, 'ref-t', fallback='NA')
            tcoupl_nvt6 = config_data.get(section, 'tcoupl', fallback='NA')
            nvt6_dir_name = f"6-NVT_re-equilibrium_integrator-{integrator_nvt6}_dt_{dt_nvt6_str}_time_{time_nvt6_ps}ps_T-{ref_t_nvt6}-K_{tcoupl_nvt6}"
            nvt6_path = os.path.join(run_parent_dir, nvt6_dir_name)
            os.makedirs(nvt6_path, exist_ok=True)
            print(f"  Criado: '{nvt6_path}'")
        else:
            print(f"Aviso: Seção '[{section}]' não encontrada. Pulando a criação deste subdiretório.")

        # --- 7. NVT Production ---
        section = '7-NVT_production'
        integrator_nvt5 = config_data.get(section, 'integrator', fallback='NA')
        dt_nvt5 = float(config_data.get(section, 'dt', fallback=0))
        dt_nvt5_str = str(dt_nvt5).replace('.', 'p')
        nsteps_nvt5 = int(config_data.get(section, 'nsteps', fallback=0))
        time_nvt5_ps = int(dt_nvt5 * nsteps_nvt5)
        ref_t_nvt5 = config_data.get(section, 'ref-t', fallback='NA')
        tcoupl_nvt5 = config_data.get(section, 'tcoupl', fallback='NA')
        nstxout_compressed_nvt5 = config_data.get(section, 'nstxout-compressed', fallback='NA')
        nvt5_dir_name = f"7-NVT_production_integrator-{integrator_nvt5}_dt_{dt_nvt5_str}_time_{time_nvt5_ps}ps_T-{ref_t_nvt5}-K_{tcoupl_nvt5}_dump-{nstxout_compressed_nvt5}"
        nvt5_path = os.path.join(run_parent_dir, nvt5_dir_name)
        os.makedirs(nvt5_path, exist_ok=True)
        print(f"  Criado: '{nvt5_path}'")
        
        return min_dir_path

    except KeyError as e:
        print(f"Erro de Chave: A variável '{e}' não foi encontrada em user_variables.txt para uma das etapas.")
        print("Por favor, verifique se todas as variáveis para os nomes dos diretórios estão definidas.")
        return None
    except ValueError as e:
        print(f"Erro de Valor: {e}. Verifique se os valores numéricos (dt, nsteps) no user_variables.txt são válidos.")
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
            # Padrão regex para encontrar o parâmetro no início da linha.
            pattern = rf"^{re.escape(param)}\s*=\s*([^;]*)(.*)$"
            
            # --- CORREÇÃO APLICADA AQUI ---
            # Adicionado re.IGNORECASE para fazer a busca sem diferenciar maiúsculas/minúsculas.
            match = re.match(pattern, line, re.IGNORECASE)
            
            if match:
                new_value = str(params_dict[param])
                comment_part = match.group(2) # O restante da linha (geralmente o comentário)

                # Mantém a formatação original da linha
                eq_pos = line.find('=')
                if eq_pos != -1:
                    param_part = line[:eq_pos].strip()
                    leading_spaces_after_param = line[len(param_part):eq_pos]
                    new_line = f"{param_part}{leading_spaces_after_param}= {new_value}{comment_part}"
                else:
                    new_line = f"{param:<25} = {new_value}{comment_part}"
                
                new_lines.append(new_line)
                replaced = True
                break
            
        if not replaced:
            new_lines.append(line)
            
    return "\n".join(new_lines)

def generate_mdp_files(main_sim_dir, config_data):
    """
    Gera os arquivos .mdp para cada etapa da simulação.
    """
    protocol_dir = os.path.join(main_sim_dir, "Protocol_input_files")
    
    md_simulation_repo_path = os.path.join(os.getcwd(), "MD_Simulation")
    template_mdp_path = os.path.join(md_simulation_repo_path, "template.mdp")

    if not os.path.exists(template_mdp_path):
        print(f"Erro: O arquivo 'template.mdp' não foi encontrado em '{template_mdp_path}'.")
        return False

    try:
        with open(template_mdp_path, 'r') as f:
            template_content = f.read()
    except IOError as e:
        print(f"Erro ao ler o arquivo 'template.mdp': {e}")
        return False

    print("\nGerando arquivos .mdp para cada etapa da simulação...")

    # Dicionário atualizado para corresponder aos novos nomes das seções
    stages = {
        "Energy Minimization": {"file_name": "1-min.mdp", "section": "1-min"},
        "NVT Thermalization":  {"file_name": "2-NVT_thermalization.mdp", "section": "2-NVT_thermalization"},
        "NPT Ergodicity":      {"file_name": "3-NPT_ergodicity.mdp", "section": "3-NPT_ergodicity"},
        "NPT Equilibration":   {"file_name": "4-NPT_equilibration.mdp", "section": "4-NPT_equilibration"},
        "NVT Re-equilibrium":  {"file_name": "6-NVT_re-equilibrium.mdp", "section": "6-NVT_re-equilibrium"}, # NOVA LINHA
        "NVT Production":      {"file_name": "7-NVT_production.mdp", "section": "7-NVT_production"},
    }


    # --- CORREÇÃO APLICADA AQUI ---
    # O nome da seção foi corrigido para 'MD_General' (com maiúsculas)
    md_general_params = {}
    if 'MD_General' in config_data:
        md_general_params = {key: value for key, value in config_data['MD_General'].items()}
    else:
        print("Aviso: Seção '[MD_General]' não encontrada em user_variables.txt.")

    success = True
    for stage_display_name, info in stages.items():
        file_name = info["file_name"]
        section_name = info["section"]
        output_mdp_path = os.path.join(protocol_dir, file_name)

        if section_name not in config_data:
            print(f"Aviso: Seção `[{section_name.upper()}]` não encontrada. Pulando a geração de '{file_name}'.")
            success = False
            continue

        stage_params_for_mdp = md_general_params.copy()
        specific_stage_config = config_data[section_name]
        
        for param, value in specific_stage_config.items():
            if value:
                stage_params_for_mdp[param] = value

        modified_content = _replace_mdp_parameters(template_content, stage_params_for_mdp)

        try:
            with open(output_mdp_path, 'w') as f:
                f.write(modified_content)
            print(f"  Gerado: '{output_mdp_path}' para '{stage_display_name}'.")
        except IOError as e:
            print(f"Erro ao escrever o arquivo '{output_mdp_path}': {e}")
            success = False
    
    return success

def _get_concentration_ratio(config_data):
    """
    Calcula a proporção das espécies moleculares em relação à primeira espécie.
    O resultado é uma string formatada como "1:2:4".

    Args:
        config_data (configparser.ConfigParser): Objeto com as variáveis do usuário.

    Returns:
        str: A string da proporção (ex: "1:2:10") ou uma string vazia se ocorrer um erro.
    """
    try:
        # Pega a quantidade da primeira espécie como linha de base para a proporção.
        # Usa getint para converter a string do arquivo para um número inteiro.
        baseline_amount = config_data.getint('general', 'specie1_ammount')

        # Se a quantidade base for zero ou não for fornecida, a proporção é indefinida.
        if baseline_amount is None or baseline_amount == 0:
            print("Erro: 'specie1_ammount' é zero ou não foi fornecido. Não é possível calcular a proporção.")
            return "RatioUndefined" # Retorna um placeholder claro

        num_species = config_data.getint('general', 'N#_species_type', fallback=0)

        ratios = []
        for i in range(1, num_species + 1):
            key = f'specie{i}_ammount'
            current_amount = config_data.getint('general', key)

            # Calcula a proporção em relação à quantidade da primeira espécie
            ratio_value = current_amount / baseline_amount

            # Formata o número para ser inteiro se não tiver casas decimais (ex: 2.0 -> 2)
            if ratio_value.is_integer():
                ratios.append(str(int(ratio_value)))
            else:
                # Se a proporção não for um inteiro, formata com 2 casas decimais
                ratios.append(f"{ratio_value:.2f}")

        return "-".join(ratios)

    except (configparser.NoOptionError, ValueError) as e:
        print(f"Erro ao ler as quantidades de espécies para calcular a proporção: {e}.")
        print("Verifique se 'N#_species_type' e todas as 'specieX_ammount' necessárias estão preenchidas com números válidos.")
        return ""
    except ZeroDivisionError:
        print("Erro: 'specie1_ammount' é zero, causando uma divisão por zero. Não é possível calcular a proporção.")
        return "RatioUndefined"
    except Exception as e:
        print(f"Um erro inesperado ocorreu em _get_concentration_ratio: {e}")
        return ""
def _get_system_size(config_data):
    """
    Soma a quantidade total de todas as espécies para definir o tamanho do sistema.
    O resultado é uma string formatada como "500species".

    Args:
        config_data (configparser.ConfigParser): Objeto com as variáveis do usuário.

    Returns:
        str: A string do tamanho total do sistema ou uma string vazia em caso de erro.
    """
    try:
        num_species = config_data.getint('general', 'N#_species_type', fallback=0)
        if num_species == 0:
            print("Aviso: 'N#_species_type' não definido em [general]. Não é possível calcular o tamanho do sistema.")
            return "0species"

        total_amount = 0
        for i in range(1, num_species + 1):
            key = f'specie{i}_ammount'
            # Usa fallback=0 para o caso de uma variável de quantidade estar vazia
            amount = config_data.getint('general', key, fallback=0)
            total_amount += amount
        
        return f"{total_amount}species"

    except (configparser.NoOptionError, ValueError) as e:
        print(f"Erro ao somar as quantidades de espécies para o tamanho do sistema: {e}.")
        print("Verifique se 'N#_species_type' e todas as 'specieX_ammount' estão preenchidas com números válidos.")
        return ""
    except Exception as e:
        print(f"Um erro inesperado ocorreu em _get_system_size: {e}")
        return ""

def _get_system_name(config_data):
    """
    Gera o nome do sistema concatenando os nomes das espécies definidos no user_variables.
    Usa 'N#_species_type' para determinar quantas espécies ler.

    Args:
        config_data (configparser.ConfigParser): Objeto com as variáveis do usuário.

    Returns:
        str: O nome do sistema formatado (ex: "c2mim-tf2n"), ou uma string vazia se ocorrer erro.
    """
    try:
        # Usamos .getint() para ler o número de forma segura
        num_species = config_data.getint('general', 'N#_species_type', fallback=0)
        if num_species == 0:
            print("Aviso: 'N#_species_type' não definido ou é zero em [general]. Não é possível gerar o nome do sistema.")
            return ""

        species_names = []
        for i in range(1, num_species + 1):
            key = f'specie{i}_name'
            # Usamos .get() para obter o nome da espécie, removendo espaços em branco
            name = config_data.get('general', key, fallback='').strip()
            if name:  # Adiciona à lista apenas se o nome não estiver vazio
                species_names.append(name)
        
        if not species_names:
            print("Aviso: Nenhuma 'specieX_name' foi encontrada, embora 'N#_species_type' seja maior que zero.")
            return ""

        # Junta os nomes com um hífen
        return "-".join(species_names)

    except (configparser.NoOptionError, ValueError) as e:
        print(f"Erro ao ler as variáveis de espécie para gerar o nome do sistema: {e}")
        return ""
    except Exception as e:
        print(f"Um erro inesperado ocorreu em _get_system_name: {e}")
        return ""

def create_subdirectories(parent_dir, config_data):
    """
    Cria os subdiretórios Force_Field_input_files, Protocol_input_files, BOX e RUN_{...}
    dentro do diretório pai.
    Retorna o caminho do diretório de minimização (ou None em caso de falha)
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
        # --- MUDANÇA AQUI ---
        return None, None, None, None, None

    print(f"\nCriando subdiretórios em '{parent_dir}'...")
    
    # ... (o resto da função permanece igual até o bloco try) ...

    fixed_subdirs = ["Force_Field_input_files", "Protocol_input_files", "BOX"]

    for subdir in fixed_subdirs:
        path = os.path.join(parent_dir, subdir)
        try:
            os.makedirs(path, exist_ok=True)
            print(f"  Criado: '{path}'")
            if subdir == "Protocol_input_files":
                protocol_path = path
        except OSError as e:
            print(f"Erro ao criar o diretório '{path}': {e}")
            # --- MUDANÇA AQUI ---
            return None, None, None, None, None

    try:
        system = _get_system_name(config_data)
        if not system:
            print("Erro: Falha ao gerar o nome do sistema. Encerrando.")
            # --- MUDANÇA AQUI ---
            return None, None, None, None, None
            
        salt_concentration = _get_concentration_ratio(config_data)
        if not salt_concentration:
            print("Aviso: Não foi possível gerar a string de concentração. Usando 'NoRatio' como placeholder.")
            salt_concentration = "NoRatio" 

        system_size = _get_system_size(config_data)
        if not system_size:
            print("Aviso: Não foi possível calcular o tamanho do sistema. Usando 'NoSize' como placeholder.")
            system_size = "NoSize"

        t_production = config_data.get('general', 'T_production', fallback='N_A') 

        run_dir_name = f"RUN_{system}_{salt_concentration}_{system_size}_{t_production}"
        run_path = os.path.join(parent_dir, run_dir_name)
        os.makedirs(run_path, exist_ok=True)
        print(f"  Criado: '{run_path}'")
        run_parent_dir = run_path

    except KeyError as e:
        print(f"Erro: Variável '{e}' necessária para o nome do diretório RUN não encontrada ou vazia.")
        # --- MUDANÇA AQUI ---
        return None, None, None, None, None
    except Exception as e:
        print(f"Um erro inesperado ocorreu ao criar o diretório RUN: {e}")
        # --- MUDANÇA AQUI ---
        return None, None, None, None, None
    
    if run_parent_dir:
        min_dir_path_return = create_run_subdirectories_actual(run_parent_dir, config_data)
        if not min_dir_path_return:
            print("Não foi possível criar todos os subdiretórios de etapa de RUN. Encerrando.")
            # --- MUDANÇA AQUI ---
            return None, None, None, None, None
    
    # --- MUDANÇA PRINCIPAL AQUI ---
    # Retornando as novas variáveis para a main()
    return min_dir_path_return, run_parent_dir, system, system_size, salt_concentration

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

    notes_dir_path = os.path.join(main_sim_dir, "User_Notes_References")

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

# Funções copy_cpu_template_job e copy_pdb_files foram removidas.

def generate_pbs_jobs(config_data, run_parent_dir, system_name, system_size, salt_concentration):
    """
    Gera os arquivos de submissão de job do PBS.
    """
    print("\nGerando arquivos de submissão de job do PBS...")

    try:
        # --- Parte 1: Ler variáveis do [PBS_job] e topology_name ---
        pbs_section = 'PBS_job'
        jobname = config_data.get(pbs_section, 'jobname', fallback='gromacs_sim')
        queue = config_data.get(pbs_section, 'queue', fallback='default')
        walltime = config_data.get(pbs_section, 'walltime', fallback='24:00:00')
        ncpus = int(config_data.get(pbs_section, 'ncpus', fallback=1))
        n_openmp = int(config_data.get(pbs_section, 'n_openmp', fallback=1))
        
        topology_name = config_data.get('general', 'topology_name', fallback='TOPOL.TOP')
        print(f"Usando topology_name: {topology_name}")


        if n_openmp > 0:
            n_mpi = ncpus // n_openmp
        else:
            n_mpi = ncpus

        # Obter nsteps e dt para 4-NPT_equilibration para o cálculo do tempo
        try:
            dt_npt4 = float(config_data.get('4-NPT_equilibration', 'dt', fallback='0'))
            nsteps_npt4 = int(config_data.get('4-NPT_equilibration', 'nsteps', fallback='0'))
            half_nsteps_ps_npt4 = (nsteps_npt4 / 2) * dt_npt4
            half_nsteps_ps_npt4_str = f"{half_nsteps_ps_npt4:.2f}".rstrip('0').rstrip('.') if '.' in str(half_nsteps_ps_npt4) else str(int(half_nsteps_ps_npt4))
            print(f"Calculado half-nsteps para 4-NPT_equilibration: {half_nsteps_ps_npt4_str} ps")
        except ValueError:
            print("Aviso: 'dt' ou 'nsteps' para 4-NPT_equilibration não são números válidos. Usando '0' para half-nsteps.")
            half_nsteps_ps_npt4_str = '0'


        # --- Parte 2: Recalcular nomes dinâmicos necessários ---
        packmol_output_pdb = f"box-{system_name}-{system_size}-species-{salt_concentration}-ratio.pdb"
        
        # Minimização 1
        emtol_em = config_data.get('1-min', 'emtol', fallback='NA')
        emtol_em_str = emtol_em.replace('.', 'p')
        emstep_em = config_data.get('1-min', 'emstep', fallback='NA')
        emstep_em_str = emstep_em.replace('.', 'p')
        min_dir_name = f"1-min_integrator-{config_data.get('1-min', 'integrator', fallback='NA')}_emtol-{emtol_em_str}_nsteps-{config_data.get('1-min', 'nsteps', fallback='NA')}_emstep-{emstep_em_str}"

        # NVT 2
        dt_nvt2 = float(config_data.get('2-NVT_thermalization', 'dt', fallback=0))
        dt_nvt2_str = str(dt_nvt2).replace('.', 'p')
        nsteps_nvt2 = int(config_data.get('2-NVT_thermalization', 'nsteps', fallback=0))
        nvt_dir_name = f"2-NVT_thermalization_integrator-{config_data.get('2-NVT_thermalization', 'integrator', fallback='NA')}_dt_{dt_nvt2_str}_time_{int(dt_nvt2 * nsteps_nvt2)}ps_T-{config_data.get('2-NVT_thermalization', 'ref-t', fallback='NA')}_{config_data.get('2-NVT_thermalization', 'tcoupl', fallback='NA')}"

        # NPT 3
        dt_npt3 = float(config_data.get('3-NPT_ergodicity', 'dt', fallback=0))
        dt_npt3_str = str(dt_npt3).replace('.', 'p')
        nsteps_npt3 = int(config_data.get('3-NPT_ergodicity', 'nsteps', fallback=0))
        ref_p_npt3 = config_data.get('3-NPT_ergodicity', 'ref-p', fallback='NA')
        ref_p_npt3_str = ref_p_npt3.replace('.', 'p')
        npt3_dir_name = f"3-NPT_ergodicity_integrator-{config_data.get('3-NPT_ergodicity', 'integrator', fallback='NA')}_dt_{dt_npt3_str}_time_{int(dt_npt3 * nsteps_npt3)}ps_T-{config_data.get('3-NPT_ergodicity', 't-initial', fallback='NA')}-{config_data.get('3-NPT_ergodicity', 't-final', fallback='NA')}-K_P-{ref_p_npt3_str}-bar_{config_data.get('3-NPT_ergodicity', 'tcoupl', fallback='NA')}_{config_data.get('3-NPT_ergodicity', 'pcoupl', fallback='NA')}"

        # NPT 4
        dt_npt4_str = str(dt_npt4).replace('.', 'p')
        ref_p_npt4 = config_data.get('4-NPT_equilibration', 'ref-p', fallback='NA')
        ref_p_npt4_str = ref_p_npt4.replace('.', 'p')
        npt4_dir_name = f"4-NPT_equilibration_integrator-{config_data.get('4-NPT_equilibration', 'integrator', fallback='NA')}_dt_{dt_npt4_str}_time_{int(dt_npt4 * nsteps_npt4)}ps_T-{config_data.get('4-NPT_equilibration', 'ref-t', fallback='NA')}-K_P-{ref_p_npt4_str}-bar_{config_data.get('4-NPT_equilibration', 'tcoupl', fallback='NA')}_{config_data.get('4-NPT_equilibration', 'pcoupl', fallback='NA')}"

        # NVT 6
        nvt6_dir_name = ""
        section_nvt6 = '6-NVT_re-equilibrium'
        if section_nvt6 in config_data:
            dt_nvt6 = float(config_data.get(section_nvt6, 'dt', fallback=0))
            dt_nvt6_str = str(dt_nvt6).replace('.', 'p')
            nsteps_nvt6 = int(config_data.get(section_nvt6, 'nsteps', fallback=0))
            nvt6_dir_name = f"6-NVT_re-equilibrium_integrator-{config_data.get(section_nvt6, 'integrator', fallback='NA')}_dt_{dt_nvt6_str}_time_{int(dt_nvt6 * nsteps_nvt6)}ps_T-{config_data.get(section_nvt6, 'ref-t', fallback='NA')}-K_{config_data.get(section_nvt6, 'tcoupl', fallback='NA')}"

        # NVT 7
        section_nvt7 = '7-NVT_production'
        dt_nvt7 = float(config_data.get(section_nvt7, 'dt', fallback=0))
        dt_nvt7_str = str(dt_nvt7).replace('.', 'p')
        nsteps_nvt7 = int(config_data.get(section_nvt7, 'nsteps', fallback=0))
        nvt7_dir_name = f"7-NVT_production_integrator-{config_data.get(section_nvt7, 'integrator', fallback='NA')}_dt_{dt_nvt7_str}_time_{int(dt_nvt7 * nsteps_nvt7)}ps_T-{config_data.get(section_nvt7, 'ref-t', fallback='NA')}-K_{config_data.get(section_nvt7, 'tcoupl', fallback='NA')}_dump-{config_data.get(section_nvt7, 'nstxout-compressed', fallback='NA')}"

        # --- Bloco 3: Gerar e Salvar cada Jobfile ---
        # --- Job 1 ---
        job1_template = f"""#!/bin/bash
#PBS -N MD1-{jobname}
#PBS -e job.md.1.err
#PBS -o job.md.1.o
#PBS -q {queue}
#PBS -l select=1:ncpus={ncpus}
#PBS -l walltime={walltime}

#EXPORT 

export I_MPI_FABRICS=shm
export OMP_NUM_THREADS={n_openmp}

#MODULES LOAD

module load intel/intel_ipp_intel64/2022.1 
module load intel/intel_ippcp_intel64/2025.1 
module load intel/mkl/2025.1 
module load intel/mpi/2021.15 
module load intel/compiler-intel-llvm/2025.1.0

#GROMASC CPU SOURCE 

source /dados/softwares/gromacs/bins/GMXRC-CPU

##############################################

cd "$PBS_O_WORKDIR"

# Job information log
{{
echo "---------------------------------------------"
echo "Job started at: $(date)"
echo "JOB ID: $PBS_JOBID"
echo "Job Name: $PBS_JOBNAME"
echo "User: $PBS_O_LOGNAME"
echo "Queue: $PBS_QUEUE"
echo "Submission Directory: $PBS_O_WORKDIR"
echo "Running Directory: $(pwd)"
echo "Running Node:"
uniq "$PBS_NODEFILE"
echo "Complete list of nodes allocated:"
cat "$PBS_NODEFILE"
echo "---------------------------------------------"
}} >> JOB.md.1.out

# Here the GROMACS part Start
# Execute GROMACS using mpirun (ensure the binary is executable)

####################1-min Running########################################################################

cd {min_dir_name}

mpirun -machinefile $PBS_NODEFILE -np 1 gmx_mpi grompp -p ../../Force_Field_input_files/{topology_name} -c ../../BOX/{packmol_output_pdb} -f ../../Protocol_input_files/1-min.mdp -o 1-min

mpirun -machinefile $PBS_NODEFILE -np {n_mpi} gmx_mpi mdrun -deffnm 1-min -ntomp {n_openmp} > mdrun.out

####################Cross-check section########################################################################

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Potential' | gmx_mpi energy -f 1-min.edr -o 1-min-potential-energy.xvg -xvg none > 1-min-potential-energy.out

tail -n 8 1-min.log > 1-min-convergence.out

####################Diretories Organization Section########################################################################

mkdir 1-min-potential-results

cp *xvg *out 1-min-potential-results/

mkdir 1-min-potential-repository

cp *edr *tpr *gro *cpt 1-min-potential-repository/

cd ..

####################2-NVT_thermalization Running########################################################################

cd {nvt_dir_name}

mpirun -machinefile $PBS_NODEFILE -np 1 gmx_mpi grompp -p ../../Force_Field_input_files/{topology_name} -f ../../Protocol_input_files/2-NVT_thermalization.mdp -c ../{min_dir_name}/1-min.gro -o 2-NVT_thermalization

mpirun -machinefile $PBS_NODEFILE -np {n_mpi} gmx_mpi mdrun -deffnm 2-NVT_thermalization -ntomp {n_openmp} > mdrun.out

####################Cross-check section######################################################################## 

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Potential' | gmx_mpi energy -f 2-NVT_thermalization.edr -o 2-NVT_thermalization-potential-energy.xvg -xvg none > 2-NVT_thermalization-potential-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Total-Energy' | gmx_mpi energy -f 2-NVT_thermalization.edr -o 2-NVT_thermalization-total-energy.xvg -xvg none > 2-NVT_thermalization-total-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Temperature' | gmx_mpi energy -f 2-NVT_thermalization.edr -o 2-NVT_thermalization-temperature.xvg -xvg none > 2-NVT_thermalization-temperature.out

####################Diretories Organization Section########################################################################

mkdir 2-NVT_thermalization-results

cp *xvg *out 2-NVT_thermalization-results/

mkdir 2-NVT_thermalization-repository

cp *edr *tpr *gro *cpt 2-NVT_thermalization-repository/

cd ..

####################3-NPT_ergodicity Running########################################################################

cd {npt3_dir_name}

mpirun -machinefile $PBS_NODEFILE -np 1 gmx_mpi grompp -p ../../Force_Field_input_files/{topology_name} -f ../../Protocol_input_files/3-NPT_ergodicity.mdp -c ../{nvt_dir_name}/2-NVT_thermalization.gro -o 3-NPT_ergodicity

mpirun -machinefile $PBS_NODEFILE -np {n_mpi} gmx_mpi mdrun -deffnm 3-NPT_ergodicity -ntomp {n_openmp} > mdrun.out

####################Cross-check section######################################################################## 

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Potential' | gmx_mpi energy -f 3-NPT_ergodicity.edr -o 3-NPT_ergodicity-potential-energy.xvg -xvg none > 3-NPT_ergodicity-potential-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Total-Energy' | gmx_mpi energy -f 3-NPT_ergodicity.edr -o 3-NPT_ergodicity-total-energy.xvg -xvg none > 3-NPT_ergodicity-total-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Temperature' | gmx_mpi energy -f 3-NPT_ergodicity.edr -o 3-NPT_ergodicity-temperature.xvg -xvg none > 3-NPT_ergodicity-temperature.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Density' | gmx_mpi energy -f 3-NPT_ergodicity.edr -o 3-NPT_ergodicity-density.xvg -xvg none > 3-NPT_ergodicity-density.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Volume' | gmx_mpi energy -f 3-NPT_ergodicity.edr -o 3-NPT_ergodicity-volume.xvg -xvg none > 3-NPT_ergodicity-volume.out

####################Diretories Organization Section########################################################################

mkdir 3-NPT_ergodicity-results

cp *xvg *out 3-NPT_ergodicity-results/

mkdir 3-NPT_ergodicity-repository

cp *edr *tpr *gro *cpt 3-NPT_ergodicity-repository/

date >> ../JOB.md.1.out
"""
        output_path1 = os.path.join(run_parent_dir, "job.md.1")
        with open(output_path1, 'w') as f:
            f.write(job1_template)
        print(f"Arquivo de job 'job.md.1' gerado com sucesso em: '{output_path1}'")

        # --- Job 2 ---
        job2_template = f"""#!/bin/bash
#PBS -N MD2-{jobname}
#PBS -e job.md.2.err
#PBS -o job.md.2.o
#PBS -q {queue}
#PBS -l select=1:ncpus={ncpus}
#PBS -l walltime={walltime}

#EXPORT 

export I_MPI_FABRICS=shm
export OMP_NUM_THREADS={n_openmp}

#MODULES LOAD

module load intel/intel_ipp_intel64/2022.1 
module load intel/intel_ippcp_intel64/2025.1 
module load intel/mkl/2025.1 
module load intel/mpi/2021.15 
module load intel/compiler-intel-llvm/2025.1.0

#GROMASC CPU SOURCE 

source /dados/softwares/gromacs/bins/GMXRC-CPU

##############################################

cd "$PBS_O_WORKDIR"

# Job information log
{{
echo "---------------------------------------------"
echo "Job started at: $(date)"
echo "JOB ID: $PBS_JOBID"
echo "Job Name: $PBS_JOBNAME"
echo "User: $PBS_O_LOGNAME"
echo "Queue: $PBS_QUEUE"
echo "Submission Directory: $PBS_O_WORKDIR"
echo "Running Directory: $(pwd)"
echo "Running Node:"
uniq "$PBS_NODEFILE"
echo "Complete list of nodes allocated:"
cat "$PBS_NODEFILE"
echo "---------------------------------------------"
}} >> JOB.md.1.out

# Here the GROMACS part Start
# Execute GROMACS using mpirun (ensure the binary is executable)

####################4-NPT_equilibration Running########################################################################

cd {npt4_dir_name}

mpirun -machinefile $PBS_NODEFILE -np 1 gmx_mpi grompp -p ../../Force_Field_input_files/{topology_name} -f ../../Protocol_input_files/4-NPT_equilibration.mdp -c ../{npt3_dir_name}/3-NPT_ergodicity.gro -o 4-NPT_equilibration

mpirun -machinefile $PBS_NODEFILE -np {n_mpi} gmx_mpi mdrun -deffnm 4-NPT_equilibration -ntomp {n_openmp} > mdrun.out

####################Cross-check section######################################################################## 

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Potential' | gmx_mpi energy -f 4-NPT_equilibration.edr -b {half_nsteps_ps_npt4_str} -o 4-NPT_equilibration-potential-energy.xvg -xvg none > 4-NPT_equilibration-potential-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Total-Energy' | gmx_mpi energy -f 4-NPT_equilibration.edr -b {half_nsteps_ps_npt4_str} -o 4-NPT_equilibration-total-energy.xvg -xvg none > 4-NPT_equilibration-total-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Temperature' | gmx_mpi energy -f 4-NPT_equilibration.edr -b {half_nsteps_ps_npt4_str} -o 4-NPT_equilibration-temperature.xvg -xvg none > 4-NPT_equilibration-temperature.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Density' | gmx_mpi energy -f 4-NPT_equilibration.edr -b {half_nsteps_ps_npt4_str} -o 4-NPT_equilibration-density.xvg -xvg none > 4-NPT_equilibration-density.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Volume' | gmx_mpi energy -f 4-NPT_equilibration.edr -b {half_nsteps_ps_npt4_str}  -o 4-NPT_equilibration-volume.xvg -xvg none > 4-NPT_equilibration-volume.out

####################Diretories Organization Section########################################################################

mkdir 4-NPT_equilibration-results

cp *xvg *out 4-NPT_equilibration-results/

mkdir 4-NPT_equilibration-repository

cp *edr *tpr *gro *cpt 4-NPT_equilibration-repository/

date >> ../JOB.md.2.out
"""
        output_path2 = os.path.join(run_parent_dir, "job.md.2")
        with open(output_path2, 'w') as f:
            f.write(job2_template)
        print(f"Arquivo de job 'job.md.2' gerado com sucesso em: '{output_path2}'")

        # --- Job 3 ---
        if nvt6_dir_name:
            job3_template = f"""#!/bin/bash
#PBS -N MD3-{jobname}
#PBS -e job.md.3.err
#PBS -o job.md.3.o
#PBS -q {queue}
#PBS -l select=1:ncpus={ncpus}
#PBS -l walltime={walltime}

#EXPORT 

export I_MPI_FABRICS=shm
export OMP_NUM_THREADS={n_openmp}

#MODULES LOAD

module load intel/intel_ipp_intel64/2022.1 
module load intel/intel_ippcp_intel64/2025.1 
module load intel/mkl/2025.1 
module load intel/mpi/2021.15 
module load intel/compiler-intel-llvm/2025.1.0

#GROMASC CPU SOURCE 

source /dados/softwares/gromacs/bins/GMXRC-CPU

##############################################

cd "$PBS_O_WORKDIR"

# Job information log
{{
echo "---------------------------------------------"
echo "Job started at: $(date)"
echo "JOB ID: $PBS_JOBID"
echo "Job Name: $PBS_JOBNAME"
echo "User: $PBS_O_LOGNAME"
echo "Queue: $PBS_QUEUE"
echo "Submission Directory: $PBS_O_WORKDIR"
echo "Running Directory: $(pwd)"
echo "Running Node:"
uniq "$PBS_NODEFILE"
echo "Complete list of nodes allocated:"
cat "$PBS_NODEFILE"
echo "---------------------------------------------"
}} >> JOB.md.1.out

# Here the GROMACS part Start
# Execute GROMACS using mpirun (ensure the binary is executable)

####################6-NVT_re-equilibrium Running########################################################################

cd {nvt6_dir_name}

mpirun -machinefile $PBS_NODEFILE -np 1 gmx_mpi grompp -p ../../Force_Field_input_files/{topology_name} -f ../../Protocol_input_files/6-NVT_re-equilibrium.mdp -c ../{npt4_dir_name}/4-NPT_equilibration.gro -o 6-NVT_re-equilibrium

mpirun -machinefile $PBS_NODEFILE -np {n_mpi} gmx_mpi mdrun -deffnm 6-NVT_re-equilibrium -ntomp {n_openmp} > mdrun.out

####################Cross-check section######################################################################## 

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Potential' | gmx_mpi energy -f 6-NVT_re-equilibrium.edr -o 6-NVT_re-equilibrium-potential-energy.xvg -xvg none > 6-NVT_re-equilibrium-potential-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Total-Energy' | gmx_mpi energy -f 6-NVT_re-equilibrium.edr -o 6-NVT_re-equilibrium-total-energy.xvg -xvg none > 6-NVT_re-equilibrium-total-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Temperature' | gmx_mpi energy -f 6-NVT_re-equilibrium.edr -o 6-NVT_re-equilibrium-temperature.xvg -xvg none > 6-NVT_re-equilibrium-temperature.out

####################Diretories Organization Section########################################################################

mkdir 6-NVT_re-equilibrium-results

cp *xvg *out 6-NVT_re-equilibrium-results/

mkdir 6-NVT_re-equilibrium-repository

cp *edr *tpr *gro *cpt 6-NVT_re-equilibrium-repository/

date >> ../JOB.md.3.out
"""
            output_path3 = os.path.join(run_parent_dir, "job.md.3")
            with open(output_path3, 'w') as f:
                f.write(job3_template)
            print(f"Arquivo de job 'job.md.3' gerado com sucesso em: '{output_path3}'")

            # --- Job 4 ---
            job4_template = f"""#!/bin/bash
#PBS -N MD4-{jobname}
#PBS -e job.md.4.err
#PBS -o job.md.4.o
#PBS -q {queue}
#PBS -l select=1:ncpus={ncpus}
#PBS -l walltime={walltime}

#EXPORT 

export I_MPI_FABRICS=shm
export OMP_NUM_THREADS={n_openmp}

#MODULES LOAD

module load intel/intel_ipp_intel64/2022.1 
module load intel/intel_ippcp_intel64/2025.1 
module load intel/mkl/2025.1 
module load intel/mpi/2021.15 
module load intel/compiler-intel-llvm/2025.1.0

#GROMASC CPU SOURCE 

source /dados/softwares/gromacs/bins/GMXRC-CPU

##############################################

cd "$PBS_O_WORKDIR"

# Job information log
{{
echo "---------------------------------------------"
echo "Job started at: $(date)"
echo "JOB ID: $PBS_JOBID"
echo "Job Name: $PBS_JOBNAME"
echo "User: $PBS_O_LOGNAME"
echo "Queue: $PBS_QUEUE"
echo "Submission Directory: $PBS_O_WORKDIR"
echo "Running Directory: $(pwd)"
echo "Running Node:"
uniq "$PBS_NODEFILE"
echo "Complete list of nodes allocated:"
cat "$PBS_NODEFILE"
echo "---------------------------------------------"
}} >> JOB.md.1.out

# Here the GROMACS part Start
# Execute GROMACS using mpirun (ensure the binary is executable)

####################7-NVT_production Running########################################################################

cd {nvt7_dir_name}

mpirun -machinefile $PBS_NODEFILE -np 1 gmx_mpi grompp -p ../../Force_Field_input_files/{topology_name} -f ../../Protocol_input_files/7-NVT_production.mdp -c ../{nvt6_dir_name}/6-NVT_re-equilibrium.gro -o 7-NVT_production

mpirun -machinefile $PBS_NODEFILE -np {n_mpi} gmx_mpi mdrun -deffnm 7-NVT_production -ntomp {n_openmp} > mdrun.out

####################Cross-check section########################################################################

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Potential' | gmx_mpi energy -f 7-NVT_production.edr -o 7-NVT_production-potential-energy.xvg -xvg none > 7-NVT_production-potential-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Total-Energy' | gmx_mpi energy -f 7-NVT_production.edr -o 7-NVT_production-total-energy.xvg -xvg none > 7-NVT_production-total-energy.out

mpirun -machinefile $PBS_NODEFILE -np 1 echo 'Temperature' | gmx_mpi energy -f 7-NVT_production.edr -o 7-NVT_production-temperature.xvg -xvg none > 7-NVT_production-temperature.out


####################Diretories Organization Section########################################################################

mkdir 7-NVT_production-results

cp *xvg *out 7-NVT_production-results/

mkdir 7-NVT_production-repository

cp *edr *tpr *gro *cpt 7-NVT_production-repository/

date >> ../JOB.md.4.out
"""
            output_path4 = os.path.join(run_parent_dir, "job.md.4")
            with open(output_path4, 'w') as f:
                f.write(job4_template)
            print(f"Arquivo de job 'job.md.4' gerado com sucesso em: '{output_path4}'")
        else:
            print("Aviso: Seção '[6-NVT_re-equilibrium]' não encontrada. Os arquivos 'job.md.3' e 'job.md.4' não serão gerados.")

        return True

    except Exception as e:
        print(f"Um erro inesperado ocorreu em generate_pbs_jobs: {e}")
        return False

def _escape_latex(text: str) -> str:
    """
    Escapa caracteres especiais do LaTeX em uma string.
    Por enquanto, focado no underscore '_'.
    """
    if not isinstance(text, str):
        text = str(text)
    return text.replace('_', '\\_')

def generate_packmol_input(config_data, system_name, system_size, salt_concentration, destination_dir):
    """
    Gera um arquivo de input do Packmol (.inp) dinamicamente com base nas
    configurações do usuário.

    Args:
        config_data (configparser.ConfigParser): Objeto com as variáveis do usuário.
        system_name (str): O nome do sistema (ex: c2mim-tf2n).
        system_size (str): O tamanho do sistema (ex: 500species).
        salt_concentration (str): A proporção das espécies (ex: 1:1).
        destination_dir (str): O diretório onde o arquivo .inp será salvo (a pasta BOX).

    Returns:
        bool: True se o arquivo foi gerado com sucesso, False caso contrário.
    """
    print("\nGerando arquivo de input do Packmol...")

    try:
        # --- Parte 1: Ler variáveis e construir os blocos das estruturas ---
        num_species = config_data.getint('general', 'N#_species_type')
        box_side = config_data.get('general', 'boxside', fallback='0.0 0.0 0.0 64.0 64.0 64.0')
        
        structure_blocks = []
        for i in range(1, num_species + 1):
            specie_name = config_data.get('general', f'specie{i}_name', fallback=None)
            specie_amount = config_data.get('general', f'specie{i}_ammount', fallback=None)

            if specie_name and specie_amount:
                # Monta o bloco para cada espécie, mantendo a indentação
                block = (
                    f"structure {specie_name}.pdb\n"
                    f"        number {specie_amount}\n"
                    f"        inside box {box_side}\n"
                    f"end structure\n"
                )
                structure_blocks.append(block)
            else:
                print(f"Aviso: 'specie{i}_name' ou 'specie{i}_ammount' não encontrados. Pulando espécie no Packmol.")

        if not structure_blocks:
            print("Erro: Nenhum bloco de estrutura foi gerado para o Packmol.")
            return False

        # --- Parte 2: Montar o conteúdo final do arquivo ---
        output_pdb_name = f"box-{system_name}-{system_size}-{salt_concentration}-ratio.pdb"
        
        header = (
            "tolerance 2.5\n"
            "filetype pdb\n"
            "randominitialpoint\n"
            "add_box_sides 0.2\n"
            f"output {output_pdb_name}\n\n"
        )

        footer = "\nseed -1\n"

        # Junta todas as partes
        final_content = header + "\n".join(structure_blocks) + footer

        # --- Parte 3: Salvar o arquivo ---
        output_filename = f"packmol-{system_name}.inp"
        output_path = os.path.join(destination_dir, output_filename)

        with open(output_path, 'w') as f:
            f.write(final_content)

        print(f"Arquivo de input do Packmol gerado com sucesso em: '{output_path}'")
        return True

    except (ValueError, configparser.NoOptionError) as e:
        print(f"Erro ao ler variáveis para o arquivo Packmol: {e}")
        return False
    except Exception as e:
        print(f"Um erro inesperado ocorreu em generate_packmol_input: {e}")
        return False


def generate_latex_log(config_data, main_sim_dir):
    """
    Gera um arquivo de log em LaTeX (log-md-protocol.tex) com uma tabela
    resumindo os principais parâmetros de cada etapa da simulação.

    Args:
        config_data (configparser.ConfigParser): Objeto com as variáveis do usuário.
        main_sim_dir (str): O caminho para o diretório principal da simulação.
    """
    print("\nGerando arquivo de log do protocolo em LaTeX...")

    # Define as seções na ordem das colunas da tabela
    sections = [
        'MD_General', '1-min', '2-NVT_thermalization', '3-NPT_ergodicity',
        '4-NPT_equilibration', '6-NVT_re-equilibrium', '7-NVT_production'
    ]
    
    # Checa quais seções realmente existem para evitar erros
    existing_sections = [s for s in sections if s in config_data or s == 'MD_General']
    
    # Mapeia cada linha da tabela para as chaves de parâmetros
    param_map = [
        ('rlist',    'integrator', 'integrator', 'integrator', 'integrator', 'integrator', 'integrator'),
        ('rcoulomb', 'emtol',      'emtol',      'emtol',      'emtol',      'emtol',      'emtol'),
        ('rvdw',     'emstep',     'emstep',     'emstep',     'emstep',     'emstep',     'emstep'),
        (None,       'tinit',      'tinit',      'tinit',      'tinit',      'tinit',      'tinit'),
        (None,       'dt',         'dt',         'dt',         'dt',         'dt',         'dt'),
        (None,       'nsteps',     'nsteps',     'nsteps',     'nsteps',     'nsteps',     'nsteps'),
        (None,       'tcoupl',     'tcoupl',     'tcoupl',     'tcoupl',     'tcoupl',     'tcoupl'),
        (None,       'tau-t',      'tau-t',      'tau-t',      'tau-t',      'tau-t',      'tau-t'),
        (None,       'ref-t',      'ref-t',      'ref-t',      'ref-t',      'ref-t',      'ref-t'),
        (None,       'gen-vel',    'gen-vel',    'gen-vel',    'gen-vel',    'gen-vel',    'gen-vel'),
        (None,       'gen-temp',   'gen-temp',   'gen-temp',   'gen-temp',   'gen-temp',   'gen-temp'),
        (None,       'pcoupl',     'pcoupl',     'pcoupl',     'pcoupl',     'pcoupl',     'pcoupl'),
        (None,       'tau-p',      'tau-p',      'tau-p',      'tau-p',      'tau-p',      'tau-p'),
        (None,       'compressibility', 'compressibility', 'compressibility', 'compressibility', 'compressibility', 'compressibility'),
        (None,       'ref-p',      'ref-p',      'ref-p',      'ref-p',      'ref-p',      'ref-p'),
        (None,       'nstlog',     'nstlog',     'nstlog',     'nstlog',     'nstlog',     'nstlog'),
        (None,       'nstenergy',  'nstenergy',  'nstenergy',  'nstenergy',  'nstenergy',  'nstenergy'),
        (None,       'nstxout-compressed', 'nstxout-compressed', 'nstxout-compressed', 'nstxout-compressed', 'nstxout-compressed', 'nstxout-compressed'),
        (None,       'xtc-precision', 'xtc-precision', 'xtc-precision', 'xtc-precision', 'xtc-precision', 'xtc-precision'),
        (None,       'continuation', 'continuation', 'continuation', 'continuation', 'continuation', 'continuation'),
        (None,       'constraints', 'constraints', 'constraints', 'constraints', 'constraints', 'constraints'),
        (None,       'constraint-algorithm', 'constraint-algorithm', 'constraint-algorithm', 'constraint-algorithm', 'constraint-algorithm', 'constraint-algorithm'),
        (None,       'T-initial',  'T-initial',  'T-initial',  'T-initial',  'T-initial',  'T-initial'),
        (None,       'T-final',    'T-final',    'T-final',    'T-final',    'T-final',    'T-final')
    ]

    fetched_values = []
    for row_keys in param_map:
        row_vals = []
        for i, key in enumerate(row_keys):
            if key and sections[i] in existing_sections:
                value = config_data.get(sections[i], key, fallback='-')
                row_vals.append(_escape_latex(value))
            else:
                row_vals.append('')
        fetched_values.append(row_vals)

    # Template LaTeX com 14 colunas agora
    latex_template = r"""
\begin{table}[]
\caption{Resumo dos parâmetros de simulação de Dinâmica Molecular utilizados em cada etapa do protocolo.}
\label{tab:md-protocol}
\resizebox{\textwidth}{!}{%
\begin{tabular}{lc|lc|lc|lc|lc|lc|lc}
\toprule
\textbf{[MD\_General]} & & \textbf{[1-min]} & & \textbf{[2-NVT\_therm.]} & & \textbf{[3-NPT\_ergod.]} & & \textbf{[4-NPT\_equilib.]} & & \textbf{[6-NVT\_re-equilib.]} & & \textbf{[7-NVT\_prod.]} & \\ \midrule
rlist & &integrator & &integrator & &integrator & &integrator & &integrator & &integrator & \\
rcoulomb & &emtol & &emtol & &emtol & &emtol & &emtol & &emtol & \\
rvdw & &emstep & &emstep & &emstep & &emstep & &emstep & &emstep & \\
 & &tinit & &tinit & &tinit & &tinit & &tinit & &tinit & \\
 & &dt & &dt & &dt & &dt & &dt & &dt & \\
 & &nsteps & &nsteps & &nsteps & &nsteps & &nsteps & &nsteps & \\
 & &tcoupl & &tcoupl & &tcoupl & &tcoupl & &tcoupl & &tcoupl & \\
 & &tau-t & &tau-t & &tau-t & &tau-t & &tau-t & &tau-t & \\
 & &ref-t & &ref-t & &ref-t & &ref-t & &ref-t & &ref-t & \\
 & &gen-vel & &gen-vel & &gen-vel & &gen-vel & &gen-vel & &gen-vel & \\
 & &gen-temp & &gen-temp & &gen-temp & &gen-temp & &gen-temp & &gen-temp & \\
 & &pcoupl & &pcoupl & &pcoupl & &pcoupl & &pcoupl & &pcoupl & \\
 & &tau-p & &tau-p & &tau-p & &tau-p & &tau-p & &tau-p & \\
 & &compressibility & &compressibility & &compressibility & &compressibility & &compressibility & &compressibility & \\
 & &ref-p & &ref-p & &ref-p & &ref-p & &ref-p & &ref-p & \\
 & &nstlog & &nstlog & &nstlog & &nstlog & &nstlog & &nstlog & \\
 & &nstenergy & &nstenergy & &nstenergy & &nstenergy & &nstenergy & &nstenergy & \\
 & &nstxout-compressed & &nstxout-compressed & &nstxout-compressed & &nstxout-compressed & &nstxout-compressed & &nstxout-compressed & \\
 & &xtc-precision & &xtc-precision & &xtc-precision & &xtc-precision & &xtc-precision & &xtc-precision & \\
 & &continuation & &continuation & &continuation & &continuation & &continuation & &continuation & \\
 & &constraints & &constraints & &constraints & &constraints & &constraints & &constraints & \\
 & &constraint-algorithm & &constraint-algorithm & &constraint-algorithm & &constraint-algorithm & &constraint-algorithm & &constraint-algorithm & \\
 & &T-initial & &T-initial & &T-initial & &T-initial & &T-initial & &T-initial & \\
 & &T-final & &T-final & &T-final & &T-final & &T-final & &T-final & \\
\bottomrule
\end{tabular}%
}
\end{table}
"""

    template_lines = latex_template.strip().splitlines()
    final_lines = []
    data_row_idx = 0
    
    for line in template_lines:
        if '&' in line and not any(cmd in line for cmd in ['\\toprule', '\\midrule', '\\bottomrule', '\\textbf']):
            if data_row_idx < len(fetched_values):
                parts = line.split('&')
                new_line = (f"{parts[0]}& {fetched_values[data_row_idx][0]} "
                            f"&{parts[2]}& {fetched_values[data_row_idx][1]} "
                            f"&{parts[4]}& {fetched_values[data_row_idx][2]} "
                            f"&{parts[6]}& {fetched_values[data_row_idx][3]} "
                            f"&{parts[8]}& {fetched_values[data_row_idx][4]} "
                            f"&{parts[10]}& {fetched_values[data_row_idx][5]} "
                            f"&{parts[12]}& {fetched_values[data_row_idx][6]} \\\\")
                final_lines.append(new_line)
                data_row_idx += 1
            else:
                 final_lines.append(line)
        else:
            final_lines.append(line)

    output_dir = os.path.join(main_sim_dir, "User_Notes_References")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "log-md-protocol.tex")

    try:
        with open(output_path, 'w') as f:
            f.write('\n'.join(final_lines))
        print(f"Arquivo de log LaTeX gerado com sucesso em: '{output_path}'")
        return True
    except IOError as e:
        print(f"Erro ao escrever o arquivo de log LaTeX '{output_path}': {e}")
        return False

def _ask_yes_no(question):
    """
    Função auxiliar que faz uma pergunta e força o usuário a responder com 'sim' ou 'não'.

    Args:
        question (str): A pergunta a ser exibida.

    Returns:
        bool: True para 'yes', False para 'não'.
    """
    while True:
        # O (s/n) no prompt indica a resposta esperada.
        response = input(f"{question} (y/n): ").lower().strip()
        if response.startswith('y'):
            return True
        elif response.startswith('n'):
            return False
        else:
            print("Invalid response. Please type 'y' for yes or 'n' for no.")

def run_interactive_checklist():
    """
    Executa um checklist interativo para o usuário no final da execução do script.
    """
    print("\n" + "="*70)
    print("--- CHECKLIST ---")
    print("="*70)
    print("System preparation is complete. Before submitting your jobs, please confirm the following points:")

    # --- Pergunta 1: Topologia e PDB ---
    question1 = (
        "\n1. Did you check your topology file and your PDB file?\n"
        "   - Is the atom order the same in both files?\n"
        "   - Are you using the correct charges for the atoms and the simulation?\n"
        "   - Remember that the simulation box MUST BE NEUTRAL.\n"
        "   - It's good practice to use the atom names in the PDB as the 'atomtypes' in the topology.\n"
    )
    if not _ask_yes_no(question1):
        print("\n!!! ATENTION: Please check your inputs again!!")
        print("The script will be terminated so you can perform the verification.")
        sys.exit()

    # --- Pergunta 2: Densidade e Volume da Caixa ---
    question2 = (
        "\n2. Did you use the experimental density to estimate the initial size of the simulation box?\n"
        "   - How did you estimate the initial box volume? Did you use PACKMOL's 'volume guesser'?\n"
    )
    if not _ask_yes_no(question2):
        print("\n!!! ATENÇÃO: Please check your inputs again!!!")
        print("The script will be terminated so you can perform the verification.")
        sys.exit()

    # --- Pergunta 3: Parâmetros do Protocolo MDP ---
    question3 = (
        "\n3. Are you 100% sure about all the parameters used in the MDP files of the protocol?\n"
        "   - Is the simulation duration appropriate for your purpose?\n"
        "   - Do the temperatures make sense for the materials?\n"
        "   - (Remember: if you are simulating a liquid, you should use a temperature between its melting and boiling points).\n"
        )
    if not _ask_yes_no(question3):
        print("\n!!! ATENÇÃO: Please check your input again!!!")
        print("O script será encerrado para que você possa fazer a verificação.")
        sys.exit()

    # --- Pergunta 4: Arquivos de Job do PBS ---
    question4 = (
        "\n4. Did you check your PBS `job.md.N` files?\n"
        "   - 'walltime', 'ncpus', '-ntomp' e '-nmpi' Are they correctly configured for the cluster queue?\n"
        "   - Remember that the product of `-nmpi` and `-ntomp` must be equal to `ncpus`."
    )
    if not _ask_yes_no(question4):
        print("\n!!! WARNING: Please check your input again!!!")
        print("The script will be terminated so you can perform the verification.")
        sys.exit()

    # --- Mensagem Final ---
    final_message = """
+----------------------------------------------------------------------+
| Good job, seems that everything is OK.                               |
|                                                                      |
| But remember to always cross-check the outputs and inputs.           |
| Always check the literature, look for similar MD works to check the  |
| protocol and workflow.                                               |
| Also update your SI_file.tex in Overleaf with the table generated    |
| here.                                                                |
|                                                                      |
| Have a nice day and good job! ✅                                     |
+----------------------------------------------------------------------+
"""
    print(final_message)
    return True

def main():
    repo_url = "git@github.com:tuananlourenco/MD_Simulation.git"
    user_variables_file = "user_variables.txt"

    if not clone_repository(repo_url):
        print("Error: Failed to clone the repository. Terminating.")
        return

    config_data = read_user_variables(user_variables_file)
    if config_data is None:
        print("Error: Failed to read the user variables file. Terminating.")
        return

    main_simulation_directory = create_main_simulation_directory(config_data)
    if not main_simulation_directory:
        print("Error: Failed to create the main simulation directory. Terminating.")
        return

    min_dir_path, run_parent_dir, system, system_size, salt_concentration = create_subdirectories(main_simulation_directory, config_data)

    if min_dir_path is None or run_parent_dir is None:
        print("Error: Failed to create essential subdirectories. Terminating.")
        return

    if not create_notes_references_directory(main_simulation_directory):
        print("Warning: Failed to create the 'User_Notes_References' directory.")

    if not generate_mdp_files(main_simulation_directory, config_data):
        print("Warning: Some MDP files could not be generated.")

    box_dir = os.path.join(main_simulation_directory, "BOX")
    if not generate_packmol_input(config_data, system, system_size, salt_concentration, box_dir):
        print("Warning: Failed to generate the Packmol input file.")

    force_field_dir = os.path.join(main_simulation_directory, "Force_Field_input_files")

    if not generate_pbs_jobs(config_data, run_parent_dir, system, system_size, salt_concentration):
        print("Warning: Failed to generate the PBS job files.")

    if not generate_latex_log(config_data, main_simulation_directory):
        print("Warning: Failed to generate the protocol log file in LaTeX format.")

    print("\nConfiguração inicial da simulação concluída. Verifique as mensagens acima para quaisquer avisos ou erros.")

    run_interactive_checklist()
    
    print("\nInitial simulation setup completed successfully.!")


if __name__ == "__main__":
    main()
