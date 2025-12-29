#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Config
# -------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="SURE-Pipe"
SURE_PIPE_SCRIPT="${SCRIPT_DIR}/SURE-Pipe"  # Launcher script
SYMLINK_PATH="/usr/local/bin/SURE-Pipe"
NEXTFLOW_DIR="$HOME/bin"

# -------------------------------
# Helper functions
# -------------------------------
command_exists() {
    command -v "$1" &>/dev/null
}
fail() { echo "❌ $1"; exit 1; }

# -------------------------------
# Conda/Mamba detection
# -------------------------------
if ! command_exists conda; then
    echo "❌ Conda is not installed. Please install Miniconda or Anaconda first."
    exit 1
fi

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

if command_exists mamba; then
    CONDA_CMD="mamba"
else
    CONDA_CMD="conda"
fi
echo "Using $CONDA_CMD for environment management."
# -------------------------------
# Install Java 17 in BASE env (conda-forge, forced)
# -------------------------------
echo "🔍 Ensuring Java 17 is correctly installed via conda-forge..."

# Add conda-forge with highest priority
conda config --add channels conda-forge 2>/dev/null || true
conda config --set channel_priority strict

echo "☕ Forcing installation of Java 17 via conda-forge..."

# Step 1: Try pulling openjdk=17
$CONDA_CMD install -y -n base -c conda-forge "openjdk>=17" --force-reinstall

# Step 2: Check Java version
JAVA_VER=$(java -version 2>&1 | awk -F[\".] '/version/ {print $2}')

if (( JAVA_VER < 17 )); then
    echo "⚠️ Java is still < 17 — resolving environment conflicts..."

    # Update all packages to allow Java 17 resolution
    $CONDA_CMD update -y -n base --all

    # Retry Java 17 installation
    $CONDA_CMD install -y -n base -c conda-forge "openjdk>=17" --force-reinstall

    JAVA_VER=$(java -version 2>&1 | awk -F[\".] '/version/ {print $2]')
fi

if (( JAVA_VER < 17 )); then
    fail "Java installation failed → Java version remains below 17."
else
    echo "☕ Java 17 successfully installed: $(java -version 2>&1 | head -n1)"
fi
# -------------------------------
# Install Nextflow if missing
# -------------------------------
if ! command_exists nextflow; then
    echo "🌱 Installing Nextflow..."
    mkdir -p "$NEXTFLOW_DIR"
    curl -s https://get.nextflow.io | bash
    mv nextflow "$NEXTFLOW_DIR/"
    chmod +x "$NEXTFLOW_DIR/nextflow"

    # Update PATH persistently
    if ! grep -q "export PATH=\"$NEXTFLOW_DIR:\$PATH\"" "$HOME/.bashrc"; then
        echo "export PATH=\"$NEXTFLOW_DIR:\$PATH\"" >> "$HOME/.bashrc"
        echo "✅ Added $NEXTFLOW_DIR to PATH in ~/.bashrc. Reload shell or run 'source ~/.bashrc'"
    fi
    export PATH="$NEXTFLOW_DIR:$PATH"
    echo "✅ Nextflow installed successfully"
else
    echo "✅ Nextflow detected: $(nextflow -v)"
fi

# -------------------------------
# Create or update Conda environment
# -------------------------------
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo "🌱 Creating ${ENV_NAME} conda environment..."
    $CONDA_CMD env create -f "${SCRIPT_DIR}/env/environment.yml"
else
    echo "♻️  Updating existing ${ENV_NAME} environment..."
    $CONDA_CMD env update -f "${SCRIPT_DIR}/env/environment.yml" --prune	
fi

# -------------------------------
# Make scripts executable
# -------------------------------
echo "🔧 Making all scripts in bin/ executable..."
if [ -d "${SCRIPT_DIR}/bin" ]; then
    find "${SCRIPT_DIR}/bin" -type f \( -name "*.sh" -o -name "*.py" \) -exec chmod +x {} \;
fi
chmod +x "$SURE_PIPE_SCRIPT"

# -------------------------------
# Symlink SURE-Pipe to PATH
# -------------------------------
if [ ! -f "$SYMLINK_PATH" ]; then
    echo "🔗 Creating symlink for SURE-Pipe in /usr/local/bin..."
    sudo ln -s "$SURE_PIPE_SCRIPT" "$SYMLINK_PATH"
else
    echo "🔗 Symlink already exists at $SYMLINK_PATH"
fi

# -------------------------------
# Install Bash Autocompletion
# -------------------------------
echo "⚙️ Installing autocompletion..."

# Ensure bash-completion is available
if ! command -v complete >/dev/null 2>&1; then
    echo "⚠️ bash-completion not found. Please install it:"
    echo "   sudo apt install bash-completion"
else
    mkdir -p ~/.bash_completion.d

    # Install completion file
    cp "${SCRIPT_DIR}/completion/surepipe.bash" ~/.bash_completion.d/surepipe.bash

    # Ensure ~/.bash_completion.d is sourced
    if ! grep -q "bash_completion.d" ~/.bashrc; then
        cat >> ~/.bashrc <<'EOF'

# Load user bash completions
if [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
fi

if [ -d "$HOME/.bash_completion.d" ]; then
    for f in "$HOME/.bash_completion.d/"*; do
        [ -f "$f" ] && source "$f"
    done
fi
EOF
    fi

    # Source immediately for current session
    if [ -f /etc/bash_completion ]; then
        source /etc/bash_completion
    fi
    for f in ~/.bash_completion.d/*; do
        [ -f "$f" ] && source "$f"
    done
fi
# -------------------------------
# Activate bashrc for current session (optional)
# -------------------------------
if [[ -n "${BASH_VERSION:-}" ]]; then
    echo "🔄 Reloading ~/.bashrc to enable autocompletion..."
    source ~/.bashrc || true

    echo "ℹ️ Bash completion installed."
    echo "➡️ Open a new terminal or run: source ~/.bashrc"
    echo "✅ Autocompletion installed and activated"
fi
