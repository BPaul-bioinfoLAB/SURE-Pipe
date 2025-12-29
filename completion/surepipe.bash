#!/usr/bin/env bash
# SURE-Pipe bash completion (final)

_surepipe_complete() {
    local cur prev
    COMPREPLY=()

    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # -------------------------------------------------
    # Locate main.nf (robust, symlink-safe)
    # -------------------------------------------------
    local main_nf=""

    # 1) Explicit override
    if [[ -n "$SUREPIPE_MAIN_NF" && -f "$SUREPIPE_MAIN_NF" ]]; then
        main_nf="$SUREPIPE_MAIN_NF"
    fi

    # 2) From installed SURE-Pipe binary
    if [[ -z "$main_nf" ]]; then
        local surepipe_bin
        surepipe_bin="$(command -v SURE-Pipe 2>/dev/null)"
        if [[ -n "$surepipe_bin" ]]; then
            surepipe_bin="$(readlink -f "$surepipe_bin")"
            local surepipe_dir
            surepipe_dir="$(cd "$(dirname "$surepipe_bin")" && pwd)"
            [[ -f "$surepipe_dir/main.nf" ]] && main_nf="$surepipe_dir/main.nf"
        fi
    fi

    # 3) Current directory fallback
    [[ -z "$main_nf" && -f "./main.nf" ]] && main_nf="./main.nf"

    [[ -f "$main_nf" ]] || return 0

    # -------------------------------------------------
    # Extract Nextflow params (all styles)
    # -------------------------------------------------
    local pipeline_params
    pipeline_params="$(
        grep -oE "params(\[['\"])?[A-Za-z0-9_]+(\['\"])?|params\.[A-Za-z0-9_]+" "$main_nf" 2>/dev/null \
        | grep -oE "[A-Za-z0-9_]+" \
        | sed 's/^/--/' \
        | sort -u
    )"

    # -------------------------------------------------
    # Static options
    # -------------------------------------------------
    local opts="run --help -h --version -V"

    local dir_params="--target_dir --neighbour_dir --output_dir --data_dir --results_dir --reference_dir"
    local file_params="--reference_genome --subject_genome --query_genome --input_genomes --genome_pairs"

    # -------------------------------------------------
    # First argument
    # -------------------------------------------------
    if [[ $COMP_CWORD -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "run ${opts}" -- "$cur") )
        return 0
    fi

    # -------------------------------------------------
    # Option name completion
    # -------------------------------------------------
    if [[ "$cur" == --* ]]; then
        COMPREPLY=( $(compgen -W "${pipeline_params} ${opts}" -- "$cur") )
        return 0
    fi

    # -------------------------------------------------
    # Directory completion
    # -------------------------------------------------
    if [[ " $dir_params " == *" $prev "* ]]; then
        compopt -o filenames 2>/dev/null
        COMPREPLY=( $(compgen -d -- "$cur") )
        return 0
    fi

    # -------------------------------------------------
    # File + directory completion
    # -------------------------------------------------
    if [[ " $file_params " == *" $prev "* || "$cur" == */* || "$prev" == */* ]]; then
        compopt -o filenames 2>/dev/null
        COMPREPLY=( $(compgen -f -- "$cur") )
        return 0
    fi

    # -------------------------------------------------
    # --param=value completion
    # -------------------------------------------------
    if [[ "$cur" == --*=* ]]; then
        local pname="${cur%%=*}"
        local pathpart="${cur#*=}"

        compopt -o filenames 2>/dev/null

        if [[ " $dir_params " == *" $pname "* ]]; then
            COMPREPLY=( $(compgen -d -- "$pathpart") )
        else
            COMPREPLY=( $(compgen -f -- "$pathpart") )
        fi

        for i in "${!COMPREPLY[@]}"; do
            COMPREPLY[$i]="${pname}=${COMPREPLY[$i]}"
        done
        return 0
    fi
}

complete -F _surepipe_complete SURE-Pipe
