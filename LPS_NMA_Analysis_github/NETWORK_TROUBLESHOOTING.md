# 🌐 GitHub 连接问题解决方案

## 🔍 问题诊断

您遇到的错误：
```
致命错误：无法访问 'https://github.com/liangjiaweihust-dot/LPS_NMA_Analysis.git/'：
Failed to connect to github.com port 443: Connection timed out
```

这通常是网络连接或防火墙问题。

## 🛠️ 解决方案

### 方案 1: 重试推送（推荐）

```bash
cd /data/lpsbinding/LPS_NMA_Analysis

# 简单重试
git push -u origin main

# 如果还是失败，尝试增加超时时间
git config --global http.timeout 300
git push -u origin main
```

### 方案 2: 使用SSH连接

```bash
# 检查是否有SSH密钥
ls -la ~/.ssh/

# 如果没有SSH密钥，生成一个
ssh-keygen -t rsa -b 4096 -C "lps.nma.analysis@example.com"

# 将SSH密钥添加到GitHub账户后，更改远程URL
git remote set-url origin git@github.com:liangjiaweihust-dot/LPS_NMA_Analysis.git

# 推送
git push -u origin main
```

### 方案 3: 配置代理（如果在企业网络中）

```bash
# 如果您知道代理服务器地址
git config --global http.proxy http://proxy.example.com:8080
git config --global https.proxy https://proxy.example.com:8080

# 推送
git push -u origin main

# 完成后可以取消代理设置
git config --global --unset http.proxy
git config --global --unset https.proxy
```

### 方案 4: 分批推送（处理大文件）

```bash
# 检查仓库大小
du -sh .git/

# 如果文件很大，尝试分批推送
git config --global http.postBuffer 524288000  # 500MB
git push -u origin main
```

### 方案 5: 使用GitHub CLI（推荐）

```bash
# 安装GitHub CLI（如果可用）
# 在CentOS/RHEL上：
sudo yum install gh

# 或者下载二进制文件
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | sudo dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null
sudo apt update
sudo apt install gh

# 认证并创建仓库
gh auth login
gh repo create LPS_NMA_Analysis --public --source=. --remote=origin --push
```

## 🔄 备用上传方案

### 方案A: 手动上传

1. **压缩项目文件**：
```bash
cd /data/lpsbinding
tar -czf LPS_NMA_Analysis.tar.gz LPS_NMA_Analysis/
```

2. **在GitHub网页端**：
   - 创建新仓库 `LPS_NMA_Analysis`
   - 选择 "uploading an existing file"
   - 拖拽压缩文件或逐个上传文件

### 方案B: 使用其他Git托管服务

```bash
# GitLab
git remote set-url origin https://gitlab.com/yourusername/LPS_NMA_Analysis.git

# Gitee (国内)
git remote set-url origin https://gitee.com/yourusername/LPS_NMA_Analysis.git

# 推送
git push -u origin main
```

### 方案C: 创建本地备份

```bash
# 创建完整备份
cd /data/lpsbinding
cp -r LPS_NMA_Analysis LPS_NMA_Analysis_backup_$(date +%Y%m%d)

# 创建发布包
cd LPS_NMA_Analysis
git archive --format=tar.gz --output=../LPS_NMA_Analysis_v1.0.0.tar.gz HEAD
```

## 🔧 网络诊断命令

```bash
# 测试GitHub连接
curl -I https://github.com

# 测试DNS解析
nslookup github.com

# 检查防火墙
sudo iptables -L | grep -i drop

# 测试端口连接
telnet github.com 443

# 检查网络路由
traceroute github.com
```

## 📋 推荐执行顺序

1. **首先尝试简单重试**：
```bash
cd /data/lpsbinding/LPS_NMA_Analysis
git config --global http.timeout 300
git push -u origin main
```

2. **如果仍然失败，检查网络**：
```bash
curl -I https://github.com
ping github.com
```

3. **考虑使用SSH**：
```bash
# 生成SSH密钥并添加到GitHub
ssh-keygen -t rsa -b 4096 -C "your.email@example.com"
cat ~/.ssh/id_rsa.pub  # 复制到GitHub SSH keys设置

# 更改为SSH URL
git remote set-url origin git@github.com:liangjiaweihust-dot/LPS_NMA_Analysis.git
git push -u origin main
```

4. **最后备选方案**：
   - 手动上传到GitHub网页
   - 使用其他Git托管服务
   - 创建本地备份

## 🎯 当前状态

您的项目已经完全准备就绪：
- ✅ Git仓库已初始化
- ✅ 所有文件已提交
- ✅ 远程仓库已配置
- ⏳ 只需要网络连接推送

**项目位置**: `/data/lpsbinding/LPS_NMA_Analysis/`  
**Git状态**: 准备推送到 `origin/main`  
**仓库大小**: 约 50-200MB（包含PDF示例文件）

## 💡 提示

如果网络问题持续存在，您可以：
1. 稍后重试（网络问题可能是临时的）
2. 联系网络管理员检查防火墙设置
3. 使用移动热点或其他网络连接
4. 考虑使用国内的Git托管服务（如Gitee）作为镜像
