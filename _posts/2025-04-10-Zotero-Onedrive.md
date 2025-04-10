---
title: Zotero + OneDrive 文献管理同步设置
date: 2025-04-10 15:16:12 +0800
categories: [Zotero, OneDrive]
tags: [Zotero, OneDrive]     # TAG names should always be lowercase
image: 
     path: /2025-04-10-Zotero-Onedrive/00.png  
---

**界面截图主要以 Mac 为例，Windows 可供参考，原理是一样的：**

1. *Zotero* -> *Setting* -> *Sync*，登录自己的账号，取消勾选 `File Syncing` 的两个选框。

     ![](/2025-04-10-Zotero-Onedrive/01.png)

2. *Zotero* -> *Setting* -> *Advanced* -> *Files and Folders* -> *Linked Attachment Base Directory*，选择 Onedrive 下的一个文件夹。（需要自己新建一个：比如 `C:\Users\XXX\OneDrive-XXX\Zotero\storage`）

     ![](/2025-04-10-Zotero-Onedrive/02.png)

3. *Data Directory Location* 的目录可以保持默认，也可以自定义个目录，建议默认。（这一步可以不用修改）

4. 打开第三步中的目录，将原有的 `sorage` 文件夹**剪切**到刚才设置的 Onedrive 文件夹中，即剪切到 `C:\Users\XXX\OneDrive-XXX\Zotero\storage` 中。（注意：第三步中的 Zotero 文件夹下不要出现 storage 文件夹。）

     ![](/2025-04-10-Zotero-Onedrive/03.png)

5. (**windows**)以管理员身份运行 **cmd** 命令行，在命令行里输入：

     ```
     mklink /J "D:\Zotero\storage" "C:\Users\***\OneDrive-XXX\Zotero\storage"
     ```

6. (**Mac**)打开 **Terminal** 终端，在命令行里输入：

     ```
     ln -s /Users/username/Library/CloudStorage/OneDrive-XXX/Zotero/storage /Users/username/Zotero/storage
     ```

     ![](/2025-04-10-Zotero-Onedrive/04.png)

以上步骤仅同步了 storage 里面的附件文件，要想同步 styels，translators，只需要将原有目录下的 styels，translators文件夹复制到 Onedrive 里对应的地方，再次运行第五或第六步的代码即可（注意修改代码里的文件夹为对应的目录）
