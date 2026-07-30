#!/bin/bash

PROJECT_DIR="/lustre/project/kennyng/damascus-earth/DaMaSCUS-NewEarth-EVAP"
cd $PROJECT_DIR

# 确保 Output 目录存在
mkdir -p Output

# 所有质量点（只提交还没提交的，或者全部）
masses="2.0 5.0 10.0"

echo "=========================================="
echo "Submitting evaporation study jobs"
echo "Output will be in: $PROJECT_DIR/Output/"
echo "Date: $(date)"
echo "=========================================="

for mass in $masses; do
    CONFIG="configs/evaporation_study/config_mass_${mass}.cfg"
    
    if [ ! -f "$CONFIG" ]; then
        echo "⚠️ Config not found: $CONFIG"
        continue
    fi
    
    echo "Submitting: mass = $mass GeV"
    
    # 使用完整的 SLURM 参数提交
    sbatch \
        --job-name="evap_${mass}" \
        --nodes=1 \
        --ntasks=4 \
        --cpus-per-task=1 \
        --time=7-00:00:00 \
        --output="Output/job-evap_${mass}-%j.out" \
        --error="Output/job-evap_${mass}-%j.err" \
        --mail-type=END,FAIL \
        --mail-user=1155234222@link.cuhk.edu.hk \
        bin/submit.sh "$CONFIG"
    
    # 等待一下，避免 QOS 限制
    sleep 2
done

echo ""
echo "=========================================="
echo "✅ All jobs submitted!"
echo "Output files will be in: $PROJECT_DIR/Output/"
echo "Check status: squeue -u $USER"
echo "=========================================="
