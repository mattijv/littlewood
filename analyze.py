import sys
import os
import re

def get_bucket_number(content):
    return re.search('bucket #([0-9]+)', content).group(1)

def assert_all_buckets_have_data(buckets, file_count):
    for i in range(1, file_count + 1):
        if str(i) not in buckets:
            print(f"Bucket {i} has no output!")
            sys.exit(1)

def get_execution_time(content):
    return float(re.search('Done in ([0-9]+\.[0-9]+) seconds', content).group(1))

def analyze_bucket(content):
    execution_time = get_execution_time(content)
    return execution_time

def get_n_longest_execution_times(buckets, n):
    bucket_list = []

    for bucket_number, execution_time in buckets.items():
        bucket_list.append([bucket_number, execution_time])
    bucket_list.sort(key = lambda bucket: bucket[1])

    return list(reversed(bucket_list))[:n]

if (len(sys.argv) < 2):
    print("Usage: analyze.py path/to/results")
    sys.exit(0)

results_dir = os.path.abspath(sys.argv[1])

print("Analyzing .out files in" + results_dir)

output_files = []
for entry in os.listdir(results_dir):
    filename, extension = os.path.splitext(entry)
    if (extension == '.out'):
        output_files.append(entry)

file_count = len(output_files)

buckets = {}

for filename in output_files:
    with open(results_dir + "/" + filename, "r") as output_file:
        content = ''.join(output_file.readlines())
        bucket_number = get_bucket_number(content)
        buckets[bucket_number] = analyze_bucket(content)

assert_all_buckets_have_data(buckets, file_count)
print("All buckets accounted for :)")

print("25 buckets that took the longest time to process:")
for bucket in get_n_longest_execution_times(buckets, 25):
    print("\t Bucket #{0: <5} took {1} seconds.".format(bucket[0], bucket[1]))
