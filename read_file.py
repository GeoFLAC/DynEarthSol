
filename = '/home/echoi2/opt/gospl_extensions/gospl_model_ext/enhanced_model.py'
try:
    with open(filename, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if 24 <= i < 100:
                print(f"{i+1}: {line.rstrip()}")
except Exception as e:
    print(e)
