#!/bin/bash

# 定义子函数
run_blast() {
    local reference_file="$1"
    local query_file="$2"
    local output_prefix="$3"
    local num_threads="$4"
    local evalue="$5"

    # 检查参数是否足够
    if [ $# -ne 5 ]; then
        echo "Error: Invalid number of arguments."
        echo "Usage: run_blast <reference_file> <query_file> <output_prefix> <num_threads> <evalue>"
        return 1
    fi

    # 删除旧的BLAST数据库文件
    if [ -f "${reference_file}.nin" ]; then rm "${reference_file}.nin"; fi
    if [ -f "${reference_file}.nsq" ]; then rm "${reference_file}.nsq"; fi
    if [ -f "${reference_file}.nhr" ]; then rm "${reference_file}.nhr"; fi
    if [ -f "${reference_file}.njs" ]; then rm "${reference_file}.njs"; fi
    if [ -f "${reference_file}.not" ]; then rm "${reference_file}.not"; fi
    if [ -f "${reference_file}.ntf" ]; then rm "${reference_file}.ntf"; fi
    if [ -f "${reference_file}.nto" ]; then rm "${reference_file}.nto"; fi
    if [ -f "${reference_file}.ndb" ]; then rm "${reference_file}.ndb"; fi

    # 创建BLAST数据库
    makeblastdb -in "$reference_file" -out "$reference_file" -dbtype nucl

    # 删除旧的输出文件
    if [ -f "${output_prefix}" ]; then rm "${output_prefix}"; fi
    if [ -f "list.txt" ]; then rm "list.txt"; fi
    if [ -f "list_uniq.txt" ]; then rm "list_uniq.txt"; fi
    if [ -f "${output_prefix}_filtered.fasta" ]; then rm "${output_prefix}_filtered.fasta"; fi

    # 运行BLAST
    blastn -query "$query_file" -db "$reference_file" -evalue "$evalue" -out "${output_prefix}" -outfmt 6 -num_threads "$num_threads"

    # 提取并去重结果
    cut -f1 "${output_prefix}" > list.txt
    sort list.txt | uniq > list_uniq.txt

    # 使用seqkit筛选结果，输入是fasta格式，所以输出结果也为fasta格式
    seqkit grep -f list_uniq.txt "$query_file" > "${output_prefix}_filtered.fasta"
    
    ########### 运行结束后，删除中间结果 ##########################
    # 删除BLAST数据库文件
    if [ -f "${reference_file}.nin" ]; then rm "${reference_file}.nin"; fi
    if [ -f "${reference_file}.nsq" ]; then rm "${reference_file}.nsq"; fi
    if [ -f "${reference_file}.nhr" ]; then rm "${reference_file}.nhr"; fi
    if [ -f "${reference_file}.njs" ]; then rm "${reference_file}.njs"; fi
    if [ -f "${reference_file}.not" ]; then rm "${reference_file}.not"; fi
    if [ -f "${reference_file}.ntf" ]; then rm "${reference_file}.ntf"; fi
    if [ -f "${reference_file}.nto" ]; then rm "${reference_file}.nto"; fi
    if [ -f "${reference_file}.ndb" ]; then rm "${reference_file}.ndb"; fi

    # 删除输出文件
    if [ -f "${output_prefix}" ]; then rm "${output_prefix}"; fi
    if [ -f "list.txt" ]; then rm "list.txt"; fi
    if [ -f "list_uniq.txt" ]; then rm "list_uniq.txt"; fi
}

# 打印帮助信息
print_help() {
    echo "Usage: $0 <reference_file> <query_file> <output_prefix> <num_threads> <evalue>"
    echo
    echo "Arguments:"
    echo "  <reference_file>  Path to the reference FASTA file."
    echo "  <query_file>      Path to the query FASTQ file."
    echo "  <output_prefix>   Prefix for the output files."
    echo "  <num_threads>     Number of threads to use for BLAST."
    echo "  <evalue>          E-value threshold for BLAST."
    echo
    echo "Example:"
    echo "  $0 reference.fa combined_subreads.fastq combined_subreads 110 1e-5"
}

# 主程序
if [ $# -ne 5 ]; then
    print_help
    exit 1
fi

# 调用子函数
run_blast "$1" "$2" "$3" "$4" "$5"





