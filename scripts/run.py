#!/usr/bin/env python
'''
this script run scaffolder with a parameter matrix
python version 3.7
Jirawat I. (nodtem66@gmail.com) Feb., 18, 2020
'''
import glob
import json
import math
import os
import subprocess
import sys
import time

import psutil


def bytes2human(n):
    # http://code.activestate.com/recipes/578019
    # >>> bytes2human(10000)
    # '9.8K'
    # >>> bytes2human(100001221)
    # '95.4M'
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.3f%s' % (value, s)
    return "%sB" % n

def run_scaffolder(input_stl, pattern, grid, coff, thickness=None, cwd=None, options=[], config=None):
    start_time = time.time()
    peak_memory = 0
    # coff = round(math.pi*coff, 5)
    if config.get('output_stl'):
        arg = [os.path.join('..', '..', 'Scaffolder.exe'), str(input_stl), str(config['output_stl']), '-n', str(pattern), '-c', str(coff), '-g', str(grid)]
    else:
        arg = [os.path.join('..', '..', 'Scaffolder.exe'), str(input_stl), '-n', str(pattern), '-c', str(coff), '-g', str(grid)]
    if thickness:
        arg.extend(['-t', str(thickness)])
    if options and type(options) is list:
        arg.extend(options)
    else:
        arg.extend(['-m', '-q', '--format', 'csv'])
    try:
        process = subprocess.Popen(arg, cwd=cwd, shell=True)
        _ps = psutil.Process(process.pid)
        while process.poll() == None:
            # from ref: https://psutil.readthedocs.io/en/latest/#psutil.Process.memory_info
            # we inqueried "Resident Set Size (RSS)" which matches “Mem Usage” column of taskmgr.exe.
            if _ps.is_running():
                if len(_ps.children()) > 0:
                    _ps = _ps.children()[0]
                peak_memory = max(peak_memory, _ps.memory_info().rss)
                time.sleep(0.3)
            else:
                break
    except:
        pass
    elapsed = (time.time() - start_time)
    mem = bytes2human(peak_memory)
    filesize = 0
    if config.get('output_stl'):
        output_stl_path = os.path.join(cwd, config['output_stl'])
        if os.path.exists(output_stl_path):
            filestat = os.stat(output_stl_path)
            filesize = bytes2human(filestat.st_size)
            os.remove(output_stl_path)
    return (elapsed, mem, filesize)

def run_batch(config, result='time_memory.txt', result_dir=None):
    input_stl = config['input_stl']
    pattern = config['pattern']
    grid_size = config['grid_size']
    coff = config['coff']
    thickness = config['thickness']
    with open(os.path.join(result_dir, result), 'w') as out:
        if type(thickness) is not list:
            thickness = [thickness]
        total = len(pattern)*len(grid_size)*len(coff)*len(thickness)
        i = 1
        for p in pattern:
            for grid in grid_size:
                for _coff in coff:
                    for t in thickness:
                        try:
                            # tracking infomation
                            print(f'{i}/{total} name={p} grid={grid} coff={_coff} thickness={t}: ', end='')
                            out.write(f'{p},{grid},{_coff*math.pi:.6},{t},')
                            # Run process
                            elapsed, mem, filesize = run_scaffolder(input_stl, p, grid, _coff, t, cwd=result_dir, options=config.get('options'), config=config)
                            # print result
                            print(f'elapsed={elapsed:.3} mem={mem} filesize={filesize}')
                            out.write(f'{elapsed:.3},{mem},{filesize}\n')
                            out.flush()
                        except Exception as e:
                            print(e)
                        i+=1

def run_list(config, result='time_memory_.txt', result_dir=None):
    input_stl = config['input_stl']
    params = config['list']
    with open(os.path.join(result_dir, result), 'w') as out:
        total = len(params)
        i = 1
        for param in params:
            if type(param) is not dict:
                print(f"Invalid {param}")
                continue
            if not param.get('pattern') or not param.get('grid_size') or not param.get('coff'):
                print("Invalid {param}: Missing \"pattern\", \"grid_size\", \"coff\"")
                continue
            thickness = 0 if not param.get('thickness') else param['thickness']
            p = param.get('pattern')
            grid = int(param.get('grid_size'))
            coff = float(param.get('coff'))
            try:
                # tracking infomation
                print(f'{i}/{total} name={p} grid={grid} coff={coff} thickness={thickness}: ', end='')
                out.write(f'{p},{grid},{coff*math.pi:.6},{thickness},')
                # Run process
                elapsed, mem, filesize = run_scaffolder(input_stl, p, grid, coff, thickness, cwd=result_dir, options=config.get('options'), config=config)
                # print result
                print(f'elapsed={elapsed:.3} mem={mem} filesize={filesize}')
                out.write(f'{elapsed:.3},{mem},{filesize}\n')
                out.flush()
            except Exception as e:
                print(e)
            i+=1
if __name__ == '__main__':
    if len(sys.argv) == 2:
        _file = sys.argv[1]
        configname = os.path.splitext(os.path.basename(_file))[0]
        result_dir = os.path.join('results', configname)
        # check path and file
        if not os.path.exists(_file):
            print(f"{_file} didn't exist")
            sys.exit(0)
        if not os.path.isfile(_file):
            print(f"{_file} is not a file")
            sys.exit(0)
        if not os.path.exists('results'):
            os.makedirs('results')
        if not os.path.exists(result_dir):
            os.makedirs(result_dir)
        else:
            _old_files = glob.glob(os.path.join(result_dir, '*.txt'))
            _old_files.extend(glob.glob(os.path.join(result_dir, '*.csv')))
            _old_files_size = len(_old_files)
            if _old_files_size > 0:
                print(f'There are {_old_files_size} result files in directory {result_dir}')
                is_removed = input('Do you want to remove? [y/n] (default no): ')
                if is_removed.lower() == 'y':
                    for f in _old_files:
                        os.remove(f)
        # read json config file
        params = {}
        with open(_file, 'r') as f:
            try:
                params = json.load(f)
            except UnicodeDecodeError as e:
                print(f"Invalid JSON file: {_file}")
                sys.exit(0)
        if not params.get('input_stl'):
            print(f"Invalid {_file}: Missing \"input_stl\" parameter")
            sys.exit(0)
        if not params.get('output'):
            print(f"Invalid {_file}: Missing \"output\" parameter")
            sys.exit(0)
        # check input stl file
        test_file = os.path.join(result_dir, params.get('input_stl'))
        if not os.path.exists(test_file):
            print(f"Not found {test_file}: \"input_stl\" parameter is relative to the folder {result_dir}")
            sys.exit(0)
        # default options
        if not params.get('thickness'):
            params['thickness'] = [0]
        mode = 'batch' if not params.get('mode') else params['mode']
        # check the options in mode
        if mode == 'batch':
            if not params.get('coff'):
                print(f"Invalid {_file}: Missing \"coff\" parameter")
                sys.exit(0)
            if not params.get('pattern'):
                print(f"Invalid {_file}: Missing \"pattern\" parameter")
                sys.exit(0)
            if not params.get('grid_size'):
                print(f"Invalid {_file}: Missing \"grid_size\" parameter")
                sys.exit(0)
            run_batch(params, params['output'], result_dir)
        elif mode == 'list':
            if not params.get('list'):
                print(f"Invalid {_file}: Missing \"list\" parameter")
                sys.exit(0)
            run_list(params, params['output'], result_dir)
        else:
            print(f"Invalid mode={mode}")
    else:
        print(f"useage: {__file__} config.json")
