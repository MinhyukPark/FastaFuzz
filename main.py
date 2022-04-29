import multiprocessing
from multiprocessing import Pool
import os
from pathlib import Path
import random
import shlex
import subprocess

import click

K = None
N = None
DEBUG = None

COMMAND_DICT = {
    "mafft": "mafft --auto INPUT_FASTA",
    "failing_executable": "./failing_executable",
}


def get_random_alphanumeric_special_sequence(length):
    return "".join(random.choice("abcdefghijklmnopqrstuvwxyz0123456789~!@#$%^&*(){}`[]?+|/=-_',.\"<>") for i in range(length))

def get_random_alphanumeric_sequence(length):
    return "".join(random.choice("abcdefghijklmnopqrstuvwxyz0123456789") for i in range(length))

def get_random_alphabet_sequence(length):
    return "".join(random.choice("abcdefghijklmnopqrstuvwxyz") for i in range(length))

def get_random_sequence_id(length):
    return "".join(random.choice("abcdefghijklmnopqrstuvwxyz") for i in range(length))

def get_random_dna_sequence(length):
    return "".join(random.choice("atcg-") for i in range(length))

def get_random_rna_sequence(length):
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
    if(not DEBUG):
        if(exit_code == 0 or exit_code == 1):
            os.remove(input_fasta)
            os.remove(output_prefix)
            os.remove(errors_prefix)
    return out,err,exit_code

def generate_single_and_multi_line_fasta(output_prefix, errors_prefix):
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"{output_prefix}/single_and_multi_line_fasta-{i}.fasta"
        current_output = f"{output_prefix}/single_and_multi_line_fasta-{i}.out"
        current_error = f"{errors_prefix}/single_and_multi_line_fasta-{i}.err"
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

def generate_varying_length_fasta(output_prefix, errors_prefix):
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"{output_prefix}/varying_length_fasta-{i}.fasta"
        current_output = f"{output_prefix}/varying_length_fasta-{i}.out"
        current_error = f"{errors_prefix}/varying_length_fasta-{i}.err"
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

def generate_duplicate_sequence_fasta(output_prefix, errors_prefix):
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"{output_prefix}/duplicate_sequence_fasta-{i}.fasta"
        current_output = f"{output_prefix}/duplicate_sequence_fasta-{i}.out"
        current_error = f"{errors_prefix}/duplicate_sequence_fasta-{i}.err"
        with open(current_filename, "w") as f:
            args_arr.append((current_filename, current_output, current_error))
            for j in range(N):
                f.write(">")
                current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                while(current_sequence_id in sequence_id_set):
                    current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                sequence_id_set.add(current_sequence_id)
                f.write(current_sequence_id + "\n")
                random_sequence = get_random_dna_sequence(random.randint(1, 500))
                f.write(random_sequence + "\n")

                while(current_sequence_id in sequence_id_set):
                    current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                sequence_id_set.add(current_sequence_id)
                f.write(current_sequence_id + "\n")
                f.write(random_sequence + "\n")
    return args_arr

def generate_duplicate_id_fasta(output_prefix, errors_prefix):
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"{output_prefix}/duplicate_id_fasta-{i}.fasta"
        current_output = f"{output_prefix}/duplicate_id_fasta-{i}.out"
        current_error = f"{errors_prefix}/duplicate_id_fasta-{i}.err"
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
                f.write(current_sequence_id + "\n")
                f.write(get_random_dna_sequence(random.randint(1, 500)) + "\n")
    return args_arr

def generate_random_alphabet_sequence_fasta(output_prefix, errors_prefix):
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"{output_prefix}/random_alphabet_sequence_fasta-{i}.fasta"
        current_output = f"{output_prefix}/random_alphabet_sequence_fasta-{i}.out"
        current_error = f"{errors_prefix}/random_alphabet_sequence_fasta-{i}.err"
        with open(current_filename, "w") as f:
            args_arr.append((current_filename, current_output, current_error))
            for j in range(N):
                f.write(">")
                current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                while(current_sequence_id in sequence_id_set):
                    current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                sequence_id_set.add(current_sequence_id)
                f.write(current_sequence_id + "\n")
                f.write(get_random_alphabet_sequence(random.randint(1, 500)) + "\n")
    return args_arr

def generate_random_alphanumeric_sequence_fasta(output_prefix, errors_prefix):
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"{output_prefix}/random_alphanumeric_sequence_fasta-{i}.fasta"
        current_output = f"{output_prefix}/random_alphanumeric_sequence_fasta-{i}.out"
        current_error = f"{errors_prefix}/random_alphanumeric_sequence_fasta-{i}.err"
        with open(current_filename, "w") as f:
            args_arr.append((current_filename, current_output, current_error))
            for j in range(N):
                f.write(">")
                current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                while(current_sequence_id in sequence_id_set):
                    current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                sequence_id_set.add(current_sequence_id)
                f.write(current_sequence_id + "\n")
                f.write(get_random_alphanumeric_sequence(random.randint(1, 500)) + "\n")
    return args_arr

def generate_random_alphanumeric_special_sequence_fasta(output_prefix, errors_prefix):
    sequence_id_set = set()
    args_arr = []
    for i in range(K):
        current_filename = f"{output_prefix}/random_alphanumeric_special_sequence_fasta-{i}.fasta"
        current_output = f"{output_prefix}/random_alphanumeric_special_sequence_fasta-{i}.out"
        current_error = f"{errors_prefix}/random_alphanumeric_special_sequence_fasta-{i}.err"
        with open(current_filename, "w") as f:
            args_arr.append((current_filename, current_output, current_error))
            for j in range(N):
                f.write(">")
                current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                while(current_sequence_id in sequence_id_set):
                    current_sequence_id = get_random_sequence_id(random.randint(1, 79))
                sequence_id_set.add(current_sequence_id)
                f.write(current_sequence_id + "\n")
                f.write(get_random_alphanumeric_special_sequence(random.randint(1, 500)) + "\n")
    return args_arr


@click.command()
@click.option("--num-processes", required=True, type=int, help="Number of processes in the main process pool")
@click.option("--output-prefix", required=True, type=click.Path(), help="Output prefix for the main log file and fuzzed fasta files")
@click.option("--errors-prefix", required=True, type=click.Path(), help="Errors prefix for the main error file")
@click.option("--num-test-cases", required=True, type=int, help="Number of test cases")
@click.option("--num-sequences", required=True, type=int, help="Number of sequences in test case")
@click.option("--debug", required=True, is_flag=True, default=False, help="Whether to keep successful input test files")
def fasta_fuzz(num_processes, output_prefix, errors_prefix, num_test_cases, num_sequences, debug):
    global K
    global N
    global DEBUG
    K = num_test_cases
    N = num_sequences
    DEBUG = debug

    pool = Pool(num_processes)
    bug_type_arr = ["single_and_multi_line", "varying_length", "duplicate_id", "duplicate_sequence", "random_alphabet_sequence", "random_alphanumeric_special"]
    args_arr = []
    args_arr.extend(generate_single_and_multi_line_fasta(output_prefix, errors_prefix))
    args_arr.extend(generate_varying_length_fasta(output_prefix, errors_prefix))
    args_arr.extend(generate_duplicate_id_fasta(output_prefix, errors_prefix))
    args_arr.extend(generate_duplicate_sequence_fasta(output_prefix, errors_prefix))
    args_arr.extend(generate_random_alphabet_sequence_fasta(output_prefix, errors_prefix))
    args_arr.extend(generate_random_alphanumeric_sequence_fasta(output_prefix, errors_prefix))
    args_arr.extend(generate_random_alphanumeric_special_sequence_fasta(output_prefix, errors_prefix))

    with open(f"{output_prefix}/main.log", "w") as f:
        with open(f"{output_prefix}/errors.log", "w") as ferr:
            f.write(f"Initialized with {K} test cases with {N} sequences each\n")
            f.write(f"Is debug mode?: {DEBUG}\n")
            f.write(f"Testing bug types: {bug_type_arr}\n")
            ferr.write(f"Collecting information about all runs...\n")
            ferr.write(f"Initialized with {K} test cases with {N} sequences each\n")
            ferr.write(f"Collecting information about error runs only...\n")
            results = pool.map(run_mafft_wrapper, args_arr)
            for i,bug_type in enumerate(bug_type_arr):
                for j in range(K):
                    results_index = i * K + j
                    exit_code = results[results_index][2]
                    if(exit_code == 0 or exit_code == 1):
                        f.write(f"SUCCESS: {bug_type} fuzzed input number {j} exited gracefully\n")
                    else:
                        f.write(f"ERROR: {bug_type} fuzzed input number {j} exited ungracefully\n")
                        ferr.write(f"ERROR: {bug_type} fuzzed input number {j} exited ungracefully\n")
            f.write(f"Logging Complete\n")
            ferr.write(f"Logging Complete\n")


if __name__ == "__main__":
    fasta_fuzz()
