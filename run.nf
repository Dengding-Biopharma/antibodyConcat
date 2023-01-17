nextflow.enable.dsl=2

process denovo {
  container 'registry-vpc.cn-hangzhou.aliyuncs.com/neobank/denovo:1.0.1'
  pod label: 'user', value: '${params.user}'
  pod imagePullSecret: 'vpc-docker'
  pod imagePullPolicy: 'Always'
  pod automountServiceAccountToken: false
  output:
    stdout

  """
  #!/usr/bin/env python3
  import subprocess
  import shutil
  import os

  for i, file_name in enumerate("${params.input}".split(",")):
    basename = os.path.basename(file_name)
    output_file = f"output/{basename}.denovo.tsv"
    process = subprocess.Popen(
      ["python3", "./denovo.py", "--model", f"/home/denovo_models/${params.modelName}", "--input", f"input/{file_name}", "--output", output_file],
      cwd="/home",
      stdout=subprocess.PIPE,
      stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    print(stdout)
    print(stderr)
  """
}

process antibodyConcat {
  container 'registry-vpc.cn-hangzhou.aliyuncs.com/neobank/antibody-concat:latest'
  pod label: 'user', value: '${params.user}'
  pod imagePullSecret: 'vpc-docker'
  pod imagePullPolicy: 'Always'
  pod automountServiceAccountToken: false
  output:
    stdout

  """
  #!/usr/bin/env python3
  import subprocess
  import shutil
  import os

  # 创建 antibodyConcat 需要的目录结构
  # -- tsv
  #   -- 1.mgf.denovo.tsv
  #   -- 2.mgf.denovo.tsv
  #   -- 3.mgf.denovo.tsv
  # -- mgf
  #   -- 1.mgf
  #   -- 2.mgf
  #   -- 3.mgf

  tsv_dir = "/home/output/tsv"
  mgf_dir = "/home/output/mgf"
  os.makedirs(tsv_dir, exist_ok=True)
  os.makedirs(mgf_dir, exist_ok=True)

  for i, file_name in enumerate("${params.input}".split(",")):
    basename = os.path.basename(file_name)

    # 把 mgf 源文件拷到 mgf 目录
    shutil.copy(f"/home/input/{file_name}", f"{mgf_dir}/{basename}")

    # 把 tsv 文件移动 tsv 目录
    shutil.copy(f"/home/output/{basename}.denovo.tsv", f"{tsv_dir}/{basename}.denovo.tsv")

  process = subprocess.Popen(
    ["python3", "run.py", "-source", "antibodyConcat", "-source_path", tsv_dir, "-spectrum_path", mgf_dir, "-score", "0.8", "-t", "2", "-kl", "5", "-ku", "8", "-predfull_path", "/install/PredFull-master/predfull.py", "-msslash_path", "/install/msSLASH-master/bin/bruteforce"],
    cwd="/home",
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
  )
  stdout, stderr = process.communicate()
  print(stdout)
  print(stderr)
  """

}

workflow {
  antibodyConcat | view
}
