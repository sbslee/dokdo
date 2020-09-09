import os
import argparse

def MakeManifest(input_path=None, output_path=None, **kwargs):
    names = []
    fnFs = []
    fnRs = []

    for r, d, f in os.walk(input_path):
        for x in f:
            if '_R1_001.fastq' in x:
                fnFs.append(f'{r}/{x}')
                names.append(x.split('_')[0])
            elif '_R2_001.fastq' in x:
                fnRs.append(f'{r}/{x}')
            else:
                pass

    with open(output_path, 'w') as f:
        headers = ['sample-id', 'forward-absolute-filepath', 
                   'reverse-absolute-filepath']
        f.write('\t'.join(headers) + '\n')

        for i in range(len(names)):
            fields = [names[i], fnFs[i], fnRs[i]]
            f.write('\t'.join(fields) + '\n')

def main():
    commands = {'MakeManifest': MakeManifest}

    parser = argparse.ArgumentParser()

    parser.add_argument('--command')
    parser.add_argument('--input-path')
    parser.add_argument('--output-path')

    args = parser.parse_args()

    commands[args.command](**vars(args))

if __name__ == '__main__':
    main()