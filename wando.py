import os
import argparse

def make_manifest(input_path=None, output_path=None, **kwargs):
    files = {}

    for r, d, f in os.walk(input_path):
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

    with open(output_path, 'w') as f:
        headers = ['sample-id', 'forward-absolute-filepath', 
                   'reverse-absolute-filepath']
        f.write('\t'.join(headers) + '\n')

        for name in sorted(files):
            fields = [name, files[name][0], files[name][1]]
            f.write('\t'.join(fields) + '\n')

def main():
    commands = {'make-manifest': make_manifest}

    parser = argparse.ArgumentParser()

    parser.add_argument('--command')
    parser.add_argument('--input-path')
    parser.add_argument('--output-path')

    args = parser.parse_args()

    commands[args.command](**vars(args))

if __name__ == '__main__':
    main()