#!/usr/bin/env python
'''
this script merge the time_memory.txt and result file (*.txt) from scaffolder.exe
python version 3.7
Jirawat I. (nodtem66@gmail.com) Feb., 19, 2020
'''
import json
import os
import sys
from collections import namedtuple
from glob import glob

Parameter = namedtuple('Parameter', ['name', 'grid', 'coff', 'thickness'])
time = {}
memory = {}

if __name__ == '__main__':
    if len(sys.argv) == 2:
        file = sys.argv[1]
        configname = os.path.splitext(os.path.basename(file))[0]
        result_dir = os.path.join('results', configname)
        if not os.path.exists(file):
            print(f"{file} didn't exist")
            sys.exit(0)
        if not os.path.isfile(file):
            print(f"{file} is not a file")
            sys.exit(0)
        params = {}
        with open(file, 'r') as f:
            try:
                params = json.load(f)
            except UnicodeDecodeError as e:
                print(f"Invalid JSON file: {file}")
                sys.exit(0)
        if not params.get('output'):
            print(f"Invalid {file}: Missing \"output\" parameter")
            sys.exit(0)
        if not params.get('merge'):
            print(f"Invalid {file}: Missing \"merge\" parameter")
            sys.exit(0)
        os.chdir(result_dir)
        output = params.get('output')
        with open(output, 'r') as fs:
            for line in fs:
                data = line.strip().split(',')
                if len(data) == 6:
                    data[2] = str(round(float(data[2]), 5))
                    data[3] = '0' if data[3] == 'None' else data[3]
                    tm = Parameter._make(data[:4])
                    time[tm] = data[4]
                    memory[tm] = data[5]
        merged_output = os.path.splitext(output)[0] + '_merged.txt'
        with open(merged_output, 'w') as out:
            results = glob(params['merge'])
            for result in results:
                with open(result, 'r') as fs:
                    data = fs.readlines()
                    if len(data) == 2:
                        data = data[1].strip()
                        d = data.split(',')
                        if len(d) >= 10:
                            coff = round(float(d[1]), 5)
                            param = Parameter(d[0], d[4], str(coff), d[3])
                            t = time.get(param)
                            m = memory.get(param)
                            if t == None or m == None:
                                print(f'{result} is invalid')
                                print('  ', end='')
                                print(param)
                            else:
                                out.write(f'{data},{t},{m}\n')
                                out.flush()
    else:
        print(f"useage: {__file__} config.json")
