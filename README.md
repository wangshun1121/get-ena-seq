# Get ENA Seq

是时候推出针对ENA上测序数据专用的下载工具了。

此工具将[ASCPsra](https://gitee.com/wangshun1121/ASCPsra)中，下载ENA数据的部分剥离出来，并利用ENA数据库的API实现数据信息自动收集整合，最终形成了这一版代码。

工具特点：

* 封装了Aspera下载命令，高效批量下载数据
* 数据传输中断可断点续传，数据传输完成自动检查
* 同时加载数据详细信息，方便数据整理与筛选
* 可以选择不直接下载数据，而是只下载数据信息，输出下载代码

## 更新信息

## 程序部署与使用

参见本项目的[Wiki](https://gitee.com/wangshun1121/get-ena-seq/wikis/)，点击下面的链接直接跳转：

* [**程序安装与环境部署**](https://gitee.com/wangshun1121/get-ena-seq/wikis/程序安装与环境部署?sort_id=3113455)

* [**从ENA下载测序数据**](https://gitee.com/wangshun1121/get-ena-seq/wikis/下载ENA测序数据?sort_id=3113459)