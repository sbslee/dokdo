import os
import argparse

def make_manifest(i_path=None, o_path=None, **kwargs):
    files = {}

    for r, d, f in os.walk(i_path):
        for x in f:
            name = x.split('_')[0]

            if '_R1_001.fastq' in x:
                if name not in files:
                    files[name] = ['', '']
                files[name][0] = f'{r}/{x}'
            elif '_R2_001.fastq' in x:
                if name not in files:
                    files[name] = ['', '']
                files[name][1] = f'{r}/{x}'
            else:
                pass

    with open(o_path, 'w') as f:
        headers = ['sample-id', 'forward-absolute-filepath', 
                   'reverse-absolute-filepath']
        f.write('\t'.join(headers) + '\n')

        for name in sorted(files):
            fields = [name, files[name][0], files[name][1]]
            f.write('\t'.join(fields) + '\n')

def merge_metadata(i_paths=None, o_path=None, **kwargs):
    metadata = []

    for i in range(len(i_paths)):
        with open(i_paths[i]) as f:
            for j, line in enumerate(f):
                fields = line.strip().split('\t')
                if i == 0 and j < 2:
                    metadata.append(fields)

                if j < 2:
                    continue

                metadata.append(fields)

    with open(o_path, 'w') as f:
        for fields in metadata:
            f.write('\t'.join(fields) + '\n')

def main():
    commands = {'make-manifest': make_manifest,
                'merge-metadata': merge_metadata}

    parser = argparse.ArgumentParser()

    parser.add_argument('--command')
    parser.add_argument('--i-path')
    parser.add_argument('--i-paths', action='append')
    parser.add_argument('--o-path')

    args = parser.parse_args()

    commands[args.command](**vars(args))

if __name__ == '__main__':
    main()