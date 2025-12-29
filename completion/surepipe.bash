#!/usr/bin/env bash
# SURE-Pipe bash completion with file/dir completion support

_surepipe_complete() {
    local cur prev opts pipeline_params main_nf
    COMPREPLY=()

    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Path to your main.nf (fallback to home location)
    main_nf="${SUREPIPE_MAIN_NF:-$HOME/Downloads/PhD/DAC_IPAC_details/SURE-Pipe/main.nf}"

    # Basic commands
    opts="run --help -h --version -V"

    # Extract params.* from Nextflow main.nf and produce --param style names
    # use uniq to avoid duplicates
    pipeline_params=$(grep -oP 'params\.\K[A-Za-z0-9_]+' "$main_nf" 2>/dev/null | sed 's/^/--/' | sort -u)

    # List of params that expect directories (complete directories)
    dir_params="--target_dir --neighbour_dir --output_dir --target_dir --data_dir --results_dir --reference_dir"

    # List of params that expect files (complete files)
    file_params="--reference_genome --subject_genome --query_genome --input_genomes --genome_pairs"

    # Helper to complete directories
    _comp_dirs() {
        # If user typed --param=/somedir<tab>
        if [[ "$cur" == *=* ]]; then
            local prefix="${cur%%=*}="
            local pathpart="${cur#*=}"
            # expand ~ if present
            local expanded
            expanded=$(eval echo "$pathpart")
            COMPREPLY=( $(compgen -d -- "$expanded") )
            # re-prefix the suggestions with --param=
            for i in "${!COMPREPLY[@]}"; do
                COMPREPLY[$i]="${prefix}$(printf '%s' "${COMPREPLY[$i]}")"
            done
        else
            # user typed: --param <tab> or after run
            COMPREPLY=( $(compgen -d -- "$cur") )
        fi
    }

    # Helper to complete files
    _comp_files() {
        if [[ "$cur" == *=* ]]; then
            local prefix="${cur%%=*}="
            local pathpart="${cur#*=}"
            local expanded
            expanded=$(eval echo "$pathpart")
            COMPREPLY=( $(compgen -f -- "$expanded") )
            for i in "${!COMPREPLY[@]}"; do
                COMPREPLY[$i]="${prefix}$(printf '%s' "${COMPREPLY[$i]}")"
            done
        else
            COMPREPLY=( $(compgen -f -- "$cur") )
        fi
    }

    # If current token begins with -- and matches param names, suggest params
    if [[ ${cur} == --* ]]; then
        COMPREPLY=( $(compgen -W "${pipeline_params} ${opts}" -- "${cur}") )
        return 0
    fi

    # If second word is run, or first word is SURE-Pipe and user typed run already
    if [[ "${prev}" == "run" || "${COMP_CWORD}" -eq 1 ]]; then
        # Suggest subcommands or params when starting the second token
        if [[ "${COMP_CWORD}" -eq 1 ]]; then
            COMPREPLY=( $(compgen -W "run ${opts}" -- "${cur}") )
            return 0
        fi
    fi

    # If previous argument is a parameter name that expects a directory
    if [[ " ${dir_params} " == *" ${prev} "* ]]; then
        _comp_dirs
        return 0
    fi

    # If previous argument is a parameter name that expects a file
    if [[ " ${file_params} " == *" ${prev} "* ]]; then
        _comp_files
        return 0
    fi

    # Also support the case where user typed --param=/path and is pressing TAB again
    # Detect parameter=value pattern in current token, and decide file/dir completion
    if [[ "$cur" == --*=* ]]; then
        # Extract param name (before '=') and normalise
        local pname="${cur%%=*}"
        if [[ " ${dir_params} " == *" ${pname} "* ]]; then
            _comp_dirs
            return 0
        elif [[ " ${file_params} " == *" ${pname} "* ]]; then
            _comp_files
            return 0
        else
            # fallback to filename completion
            COMPREPLY=( $(compgen -f -- "${cur#*=}") )
            # prefix back
            for i in "${!COMPREPLY[@]}"; do
                COMPREPLY[$i]="${pname}=${COMPREPLY[$i]}"
            done
            return 0
        fi
    fi

    return 0
}

complete -F _surepipe_complete SURE-Pipe
