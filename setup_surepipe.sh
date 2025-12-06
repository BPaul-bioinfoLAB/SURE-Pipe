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
# Check Java version (must be >=17)
# -------------------------------
if command_exists java; then
    JAVA_VER=$(java -version 2>&1 | awk -F[\".] '/version/ {print $2}')
    if (( JAVA_VER < 17 )); then
        fail "Java version < 17 detected. Install Java 17+ using SDKMAN:
  sdk install java 17.0.10-tem"
    fi
    echo "☕ Java >=17 detected: $(java -version 2>&1 | head -n1)"
else
    fail "Java not installed. Install Java 17+ using:
  curl -s https://get.sdkman.io | bash
  sdk install java 17.0.10-tem"
fi


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
# Done
# -------------------------------
echo "✅ Environment setup complete!"
echo "Activate it with: conda activate ${ENV_NAME}"
echo "Run SURE-Pipe anywhere: SURE-Pipe <command>"

