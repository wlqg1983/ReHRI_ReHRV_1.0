#!/bin/bash

# Define color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default environment name
DEFAULT_ENV_NAME="ReHRI_ReHRV_1.0"
YML_FILE="ReHRI_ReHRV_1.0.yml"
MAX_RETRIES=3

# Function to validate environment name
validate_env_name() {
    local env_name=$1
    if [[ ! $env_name =~ ^[a-zA-Z][a-zA-Z0-9_]*$ ]]; then
        echo -e "${RED}Error: Environment name '$env_name' is invalid!${NC}"
        echo -e "${YELLOW}Environment name must:"
        echo -e "1. Start with a letter"
        echo -e "2. Contain only letters, numbers, or underscores"
        echo -e "3. Not contain spaces or special characters${NC}"
        exit 1
    fi
}

# Function to check conda installation
check_conda() {
    if ! command -v conda &> /dev/null; then
        echo -e "${RED}Error: Conda is not installed or not in PATH!${NC}"
        echo -e "${YELLOW}Please install conda first:"
        echo -e "1. Download Miniconda: https://docs.conda.io/en/latest/miniconda.html"
        echo -e "2. Install it and add to PATH"
        echo -e "3. Restart your terminal before running this script again${NC}"
        exit 1
    fi
}

# Function to attempt installation with retries
install_with_retry() {
    local install_cmd=$1
    local component=$2
    local success=0
    
    for ((i=1; i<=$MAX_RETRIES; i++)); do
        echo -e "${YELLOW}Attempt $i/$MAX_RETRIES to install $component...${NC}"
        
        if eval "$install_cmd"; then
            echo -e "${GREEN}Successfully installed $component on attempt $i${NC}"
            success=1
            break
        else
            echo -e "${YELLOW}Attempt $i failed. Retrying in 3 seconds...${NC}"
            sleep 3
        fi
    done
    
    if [ $success -eq 0 ]; then
        echo -e "${RED}Failed to install $component after $MAX_RETRIES attempts${NC}"
        return 1
    fi
    return 0
}

# Function to show manual installation instructions
show_install_help() {
    local component=$1
    local install_cmd=$2
    
    echo -e "\n${BLUE}>>>> Manual Installation Required <<<<${NC}"
    echo -e "${YELLOW}Component:${NC} $component"
    echo -e "${YELLOW}Possible issue:${NC} Network connectivity or package repository problems"
    echo -e "${YELLOW}Solution:${NC}"
    echo -e "1. Check your internet connection"
    echo -e "2. Try running the installation command manually:"
    echo -e "   ${GREEN}$install_cmd${NC}"
    echo -e "3. If using a corporate network, you may need to configure proxies"
    echo -e "4. For conda issues, try: ${GREEN}conda clean -i && conda update -n base conda${NC}"
    echo -e "5. For pip issues, try: ${GREEN}pip install --upgrade pip${NC}"
    echo -e "${BLUE}>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<${NC}\n"
}

# Handle command line arguments
if [ $# -eq 0 ]; then
    ENV_NAME=$DEFAULT_ENV_NAME
elif [ $# -eq 1 ]; then
    ENV_NAME=$1
    validate_env_name "$ENV_NAME"
else
    echo -e "${RED}Error: Invalid number of arguments!${NC}"
    echo -e "${YELLOW}Usage:"
    echo -e "  $0 [environment_name]"
    echo -e "  If no environment name provided, defaults to '$DEFAULT_ENV_NAME'${NC}"
    exit 1
fi

# Check conda installation
check_conda

# 1. Create conda environment with retries
echo -e "${YELLOW}=== Creating conda environment '$ENV_NAME'... ===${NC}"
if ! install_with_retry "conda env create -n $ENV_NAME -f $YML_FILE" "conda environment"; then
    show_install_help "Conda Environment" "conda env create -n $ENV_NAME -f $YML_FILE"
    exit 1
fi

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# 2. Install plasmidrender with retries
echo -e "\n${YELLOW}=== Installing plasmidrender... ===${NC}"
plasmidrender_path="$(conda info --base)/envs/$ENV_NAME/lib/python3.12/site-packages/plasmidrender"
if [ -d "bin/plasmidrender" ]; then
    for ((i=1; i<=$MAX_RETRIES; i++)); do
        echo -e "${YELLOW}Attempt $i/$MAX_RETRIES to copy plasmidrender...${NC}"
        if cp -r bin/plasmidrender $plasmidrender_path; then
            echo -e "${GREEN}Successfully copied plasmidrender on attempt $i${NC}"
            rm -rf bin/plasmidrender
            chmod +x bin/*
            break
        else
            echo -e "${YELLOW}Attempt $i failed. Retrying in 3 seconds...${NC}"
            sleep 3
            if [ $i -eq $MAX_RETRIES ]; then
                echo -e "${RED}Failed to copy plasmidrender after $MAX_RETRIES attempts${NC}"
                show_install_help "plasmidrender" "Manual copy required: cp -r /path/to/plasmidrender $plasmidrender_path"
            fi
        fi
    done
else
    echo -e "${RED}plasmidrender source directory not found in bin/!${NC}"
    show_install_help "plasmidrender" "Manual copy required: cp -r /path/to/plasmidrender $plasmidrender_path"
fi

# 3. Verification process with retries
echo -e "\n${YELLOW}=== Verifying installation... ===${NC}"

# Tools to verify (key=command, value=package)
declare -A TOOLS=(
    ["blastn"]="blast"
    ["samtools"]="samtools"
    ["minimap2"]="minimap2"
    ["seqkit"]="seqkit"
    ["python"]="python"
)

# Python modules to verify (key=module, value=package)
declare -A MODULES=(
    ["biopython"]="biopython"
    ["pandas"]="pandas"
    ["matplotlib"]="matplotlib"
    ["numpy"]="numpy"
    ["plasmidrender"]="plasmidrender"
)

# Track installation failures
declare -a FAILED_INSTALLS=()

# Verify command line tools
echo -e "\n${YELLOW}--- Command Line Tools ---${NC}"
for tool in "${!TOOLS[@]}"; do
    if which $tool >/dev/null; then
        echo -e "${GREEN}[✓] $tool found: $(which $tool)${NC}"
    else
        echo -e "${RED}[✗] $tool missing!${NC}"
        if ! install_with_retry "conda install -y ${TOOLS[$tool]}" "$tool"; then
            FAILED_INSTALLS+=("$tool")
            show_install_help "$tool" "conda install -y ${TOOLS[$tool]}"
        fi
    fi
done

# Verify Python modules
echo -e "\n${YELLOW}--- Python Modules ---${NC}"
for module in "${!MODULES[@]}"; do
    if python -c "import $module" 2>/dev/null; then
        echo -e "${GREEN}[✓] $module imports successfully${NC}"
    else
        echo -e "${RED}[✗] Failed to import $module!${NC}"
        
        if [[ $module == "plasmidrender" ]]; then
            if [ -d "$plasmidrender_path" ]; then
                echo -e "${GREEN}[✓] plasmidrender directory exists${NC}"
            else
                echo -e "${RED}[✗] plasmidrender installation incomplete!${NC}"
                FAILED_INSTALLS+=("$module")
                show_install_help "$module" "Manual copy required: cp -r /path/to/plasmidrender $plasmidrender_path"
            fi
        else
            if ! install_with_retry "pip install ${MODULES[$module]}" "$module"; then
                FAILED_INSTALLS+=("$module")
                show_install_help "$module" "pip install ${MODULES[$module]}"
            fi
        fi
    fi
done

# Verify bin directory
echo -e "\n${YELLOW}--- Executable Files ---${NC}"
bin_executables=$(find bin/ -type f -exec test -x {} \; -print | wc -l)
if [ $bin_executables -gt 0 ]; then
    echo -e "${GREEN}[✓] Found $bin_executables executable files${NC}"
else
    echo -e "${RED}[✗] No executable files found!${NC}"
    for ((i=1; i<=$MAX_RETRIES; i++)); do
        echo -e "${YELLOW}Attempt $i/$MAX_RETRIES to fix permissions...${NC}"
        chmod +x bin/*
        bin_executables=$(find bin/ -type f -exec test -x {} \; -print | wc -l)
        if [ $bin_executables -gt 0 ]; then
            echo -e "${GREEN}[✓] Permissions fixed successfully on attempt $i${NC}"
            break
        else
            echo -e "${YELLOW}Attempt $i failed. Retrying in 3 seconds...${NC}"
            sleep 3
            if [ $i -eq $MAX_RETRIES ]; then
                echo -e "${RED}Failed to fix permissions after $MAX_RETRIES attempts${NC}"
                FAILED_INSTALLS+=("bin_executables")
                show_install_help "Bin Directory" "chmod +x bin/*"
            fi
        fi
    done
fi

# Final report
echo -e "\n${YELLOW}=== Installation Summary ===${NC}"
if [ ${#FAILED_INSTALLS[@]} -eq 0 ]; then
    echo -e "${GREEN}All components installed successfully!${NC}"
else
    echo -e "${RED}The following components failed to install after $MAX_RETRIES attempts:${NC}"
    for item in "${FAILED_INSTALLS[@]}"; do
        echo -e " - ${RED}$item${NC}"
    done
    echo -e "\n${YELLOW}Refer to the instructions above each failure for manual installation steps."
    echo -e "Common solutions include:"
    echo -e "1. Checking your internet connection"
    echo -e "2. Running the suggested install commands manually"
    echo -e "3. Configuring proxies if behind a corporate firewall"
    echo -e "4. Updating conda/pip:"
    echo -e "   conda: ${GREEN}conda clean -i && conda update -n base conda${NC}"
    echo -e "   pip:   ${GREEN}pip install --upgrade pip${NC}${NC}"
fi

# Clean up
rm -f $YML_FILE

