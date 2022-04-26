import random
import subprocess
import multiprocessing
from pathlib import Path
from multiprocessing import Pool
import shlex

K = 30
N = 30
COMMAND_DICT = {
    "mafft": "mafft --auto INPUT_FASTA",
    "failing_executable": "./failing_executable",
}


def get_random_sequence_id(length):
    return "".join(random.choice("abcdefghijklmnopqrstuvwxyz") for i in range(length))

def get_random_dna_sequence(length):
    return "".join(random.choice("atcg-") for i in range(length))

def get_random_dna_sequence(length):
    return "".join(random.choice("augc-") for i in range(length))

def run_mafft_wrapper(args):
    return run("mafft", *args)
    # return run("failing_executable", *args)

def run(command, input_fasta, output_prefix, errors_prefix):
    command = COMMAND_DICT[command].replace("INPUT_FASTA", input_fasta)
    with open(output_prefix, "w") as fout:
        with open(errors_prefix, "w") as ferr:
            # print(shlex.split(iqtree_command))
            p = subprocess.Popen(shlex.split(command), stdout=fout, stderr=ferr)
    out,err = p.communicate()
    exit_code = p.wait()
    return out,err,exit_code

def generate_single_and_multi_line_fasta():
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"./output/single_and_multi_line_fasta-{i}.fasta"
        current_output = f"./output/single_and_multi_line_fasta-{i}.out"
        current_error = f"./errors/single_and_multi_line_fasta-{i}.err"
        with open(current_filename, "w") as f:
            args_arr.append((current_filename, current_output, current_error))
            for j in range(N):
                f.write(">")
                current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                while(current_sequence_id in sequence_id_set):
                    current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                sequence_id_set.add(current_sequence_id)
                f.write(current_sequence_id + "\n")
                if(i < N / 2):
                    f.write(get_random_dna_sequence(79) + "\n")
                else:
                    f.write(get_random_dna_sequence(79) + "\n")
                    f.write(get_random_dna_sequence(79) + "\n")
                    f.write(get_random_dna_sequence(79) + "\n")
                    f.write(get_random_dna_sequence(79) + "\n")
    return args_arr

def generate_varying_length_fasta():
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"./output/varying_length_fasta-{i}.fasta"
        current_output = f"./output/varying_length_fasta-{i}.out"
        current_error = f"./errors/varying_length_fasta-{i}.err"
        with open(current_filename, "w") as f:
            args_arr.append((current_filename, current_output, current_error))
            for j in range(N):
                f.write(">")
                current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                while(current_sequence_id in sequence_id_set):
                    current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                sequence_id_set.add(current_sequence_id)
                f.write(current_sequence_id + "\n")
                f.write(get_random_dna_sequence(random.randint(1, 500)) + "\n")
    return args_arr

def main():
    pool = Pool(16)
    args_arr = []
    args_arr.extend(generate_single_and_multi_line_fasta())
    args_arr.extend(generate_varying_length_fasta())
    # args_arr.extend(generate_duplicate_id_fasta())
    # args_arr.extend(generate_duplicate_sequence_fasta())

    results = pool.map(run_mafft_wrapper, args_arr)
    print(results)
    # results = pool.map(run_iqtree_wrapper, args_arr)
    # print(results)
    # results = pool.map(run_iqtree_wrapper, args_arr)
    # print(results)


if __name__ == "__main__":
    main()
