---
title: Ruby + Jekyll + Bundle 配置网页（Github Pages）
date: 2025-04-10 15:16:12 +0800
categories: [Ruby, Jekyll]
tags: [Ruby, Jekyll, Github]     # TAG names should always be lowercase
image: 
     path: /2025/04-11-00.png  
description: 不是经验分享。是给自己的笔记，只有自己能看懂
---

> 不是经验分享。是给自己的笔记，**只有自己能看懂**


**仅展示 Windows 环境**，Mac 部分步骤类似，命令也是相同的：


1. 初次使用，先安装 **Ruby+Devkit**： [https://rubyinstaller.org/downloads/](https://rubyinstaller.org/downloads/)
2. 安装完之后选择：```MSYS2 and MINGW development toolchain```
3. 打开新终端：```gem install jekyll bundler```
4. 测试jekyll是否安装成功：```jekyll -v```

5. 然后在blog目录下，运行 bundle，安装所需依赖；

6. 构建与部署命令：

    ```bundle exec jekyll s```

    部署到 [http://127.0.0.1:4000](http://127.0.0.1:4000)

7. 然后通过git提交：

    ```bash
    git add .
    git commit -m "updates"
    git push
    ```