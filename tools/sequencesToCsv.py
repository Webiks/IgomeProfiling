import argparse


def sequences_to_csv(input_path, output_path):
    with open(input_path, 'r') as f:
        lines = f.readlines()
    
    gaps = []
    count = 0

    with open(output_path, 'w') as f:
        for i in range(1, len(lines), 2):
            count += 1
            chars = list(lines[i].rstrip())
            if len(gaps) == 0:
                gaps = [0 for c in chars]
            for j in range(len(chars)):
                if chars[j] == '-':
                    gaps[j] += 1
            f.write(','.join(chars))
            f.write('\n')
    
        for i in range(len(gaps)):
            gaps[i] /= count

        f.write('\n')
        f.write(','.join([str(gap) for gap in gaps]))
        f.write('\n')

    print(gaps)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('faa', help='A path to faa file with sequence')
    parser.add_argument('csv', help='Output path of sequences separated by letters')
    args = parser.parse_args()
    
    sequences_to_csv(args.faa, args.csv)
