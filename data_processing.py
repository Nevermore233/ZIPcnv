from collections import defaultdict
import pandas as pd
import time
from tqdm import trange
import os
from datetime import datetime
from utils import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-config', type=str, help="Path to the '.config' file ", required=True)
args = parser.parse_args()


def standardize_depth(sample_depths, chr_len_list, epsilon=1e-6, save_filename="rj_means_and_n.json"):
    standardized_depths = defaultdict(dict)
    rj_means_dict = {}
    n_samples = len(sample_depths)

    for chr_name, chr_len in chr_len_list:
        # Collect depths for all samples at this chromosome
        all_samples_depth = np.array([sample_depth[chr_name] for sample_depth in sample_depths])

        # Compute Rj^mode for each sample, adding epsilon to avoid division by zero
        Rj_median = np.median(all_samples_depth, axis=1) + epsilon

        # Compute Ri·^mean for each position, adding epsilon to avoid division by zero
        Ri_means = np.mean(all_samples_depth, axis=0) + epsilon

        rj_means_dict[chr_name] = {
            "Rj_means": Ri_means.tolist(),  # 转换为列表格式
            "n_samples": n_samples
        }

        for sample_idx, sample_depth in enumerate(all_samples_depth):
            # Standardize depth
            standardized_depth = (sample_depth / Rj_median[sample_idx]) / Ri_means
            standardized_depths[sample_idx][chr_name] = np.round(standardized_depth, 2)

    with open(save_filename, "w") as f:
        json.dump(rj_means_dict, f, indent=4)

    return standardized_depths


def get_std_dep(df, chr_len_list, read_len, log_filename):
    sample_depths = []
    for i in trange(df.shape[0]):
        mapping = f"sample_{i}"
        filename = df.loc[df['mapping'] == mapping]['file_name'].values[0]
        print(f'process {filename} ..................')
        if os.path.isfile(filename):
            try:
                sample_depth = calcu_bam_dep(chr_len_list, filename, read_len)
                sample_depths.append(sample_depth)
            except Exception as e:
                current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
                error_message = f"{current_time}:Error occurred while reading {filename}: {e}"
                print(error_message)
                # Write the error message to the log file
                with open(log_filename, "a") as log_file:
                    log_file.write(error_message + "\n")
                continue
    return sample_depths


def save2json(standardized_depths, filename):
    # Convert numpy arrays to lists for JSON serialization
    standardized_depths_serializable = {}
    for chr_name, standardized_depth in standardized_depths.items():
        standardized_depths_serializable[chr_name] = standardized_depth.tolist()

    # Write to JSON file
    with open(filename, 'w') as jsonfile:
        json.dump(standardized_depths_serializable, jsonfile)


def main():
    # .config file
    config_file = args.config
    paths = read_config(config_file)

    # Input tested files
    test_file_list = paths['test_file_list']
    test_df = pd.read_csv(test_file_list, index_col=0)

    # Baseline files, at least 50 normal samples
    baseline_file_list = paths['baseline_file_list']
    bl_df = pd.read_csv(baseline_file_list, index_col=0)

    # Train (if exists)
    try:
        train_file_list = paths['train_file_list']
        if train_file_list == '0' or train_file_list is None:
            print("Notice! Training samples do not exist.")
            train_df = None
        else:
            train_df = pd.read_csv(train_file_list, index_col=0)
    except (FileNotFoundError, pd.errors.EmptyDataError, KeyError, ValueError) as e:
        print("Warning! Training samples do not exist.")
        train_df = None

    all_data = pd.concat([test_df['file_name'], bl_df['file_name']], ignore_index=True)

    if train_df is not None:
        all_data = pd.concat([all_data, train_df['file_name']], ignore_index=True)

    all_data = all_data.to_frame(name='file_name')
    all_data['mapping'] = ['sample_' + str(i) for i in range(len(all_data))]

    # Chromosome bed file
    chr_len_path = paths['chr_len_path']
    chr_len_list = read_chr_len_file(chr_len_path)
    read_len = int(paths['read_len'])

    # Logs
    log_filename = "log/data_processing_log.txt"
    os.makedirs(os.path.dirname(log_filename), exist_ok=True)

    # Standardization
    all_sample_depths = get_std_dep(all_data, chr_len_list, read_len, log_filename)
    all_standardized_depths_list = standardize_depth(all_sample_depths, chr_len_list)

    # Save
    print('save to file ...')

    # Create output directory if it doesn't exist
    json_dir = 'data/nor/'

    if not os.path.exists(json_dir):
        os.makedirs(json_dir)

    for i in trange(all_data.shape[0]):
        # Mapping
        mapping = f"sample_{i}"
        sample = all_data.loc[all_data['mapping'] == mapping]['file_name'].values[0]

        filename_part = os.path.basename(sample).replace('.bam', '.json')

        new_filename = os.path.join(json_dir, filename_part)

        standardized_depths = all_standardized_depths_list[i]
        save2json(standardized_depths, new_filename)

    # Set a baseline for comparison
    # Baseline save path
    baseline_save_path = paths['baseline_save_path']

    if not os.path.exists(baseline_save_path):
        os.makedirs(baseline_save_path)

    os.makedirs(os.path.dirname(log_filename), exist_ok=True)

    # Warning
    if bl_df.shape[0] < 50:
        print('WARNING: Please input at least 50 samples as a baseline.')

    cb_data = {}
    for chr_name, chr_len in chr_len_list:
        cb_data[chr_name] = np.zeros(chr_len)
    for i in trange(bl_df.shape[0]):
        mapping = f"sample_{i}"
        sample = bl_df.loc[bl_df['mapping'] == mapping]['file_name'].values[0]
        filename_part_ = os.path.basename(sample).replace('.bam', '.json')
        filename = os.path.join(json_dir, filename_part_)
        print(f'Process {filename} ..................')
        if os.path.isfile(filename):
            try:
                sample_depth = load_from_json(filename)
                for chr_name, chr_len in chr_len_list:
                    cb_data[chr_name] = cb_data[chr_name] + sample_depth[chr_name]
            except Exception as e:
                current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
                error_message = f"{current_time}:Error occurred while reading {filename}: {e}"
                print(error_message)
                # Write the error message to the log file
                with open(log_filename, "a") as log_file:
                    log_file.write(error_message + "\n")
                continue

    # Save to file
    for chr_name in cb_data:
        cb_data[chr_name] /= bl_df.shape[0]
        np.savez(f'{baseline_save_path}/baseline_file_{chr_name}', **cb_data)

    if bl_df.shape[0] < 50:
        print('WARNING: Please input at least 50 samples as a baseline.')


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"Data process finish! runtime: {rt}sec")


