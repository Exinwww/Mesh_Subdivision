# Mesh_Subdivision
Loop Subdivision by python(CG_work)


MeshSub.py 中 点存为Point类， 面为Face类；

使用Load_OBJ读入.obj文件
目前仅处理v，vn， f数据，vt并未有相应操作；

修改以下内容，运行程序：
rootPath为模型文件根目录
filePath为模型文件的相对地址
二者传入执行函数Subdivision（）中开始执行；
还可以额外传入细分次数times，指定细分运行次数；
[tips：模型文件较大时，运行次数过多会耗时过长，内存消耗大]

仍有不足，尽请指教
