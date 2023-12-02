import numpy as np
np.set_printoptions(precision=6)
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

class Point:
    def __init__(self, pos, order = None):
        self.pos = np.array(pos)
        self.order = order
    def __str__(self)->str:
        res = ['v', f'{self.pos[0]:.6f}', f'{self.pos[1]:.6f}', f'{self.pos[2]:.6f}']
        return ' '.join(res)
    def __lt__(self, other):
        assert isinstance(other, Point)
        return self.order < other.order

class Face:
    def __init__(self, v, vt = None, vn = None):
        # v is index of vertex
        v = np.array([int(ii) for ii in v])
        self.v = v
        self.vt = vt
        if vn is None: self.vn = vn
        else: self.vn = self.v
    def __str__(self)->str:
        res = ['f']
        datas = [self.v, self.vt, self.vn]

        for i in range(3):
            tString = ['', '', '']
            for ii in range(3):
                temp = datas[ii]
                if temp is None or len(temp)==0:
                    continue
                tString[ii] =  str(temp[i])
            temp = '/'.join(tString)
            if temp.endswith('//'): temp = temp[:len(temp)-2]
            res.append(temp)
        return ' '.join(res)
    def getNormal(self, vs)->np.ndarray:
        assert isinstance(vs, list)
        assert len(self.v) == 3
        # a, b, c = self.v
        # assert 1==0 #!!!!!!!!!!!!error!!!!!!!!!!
        a = vs[int(self.v[0])-1].pos
        b = vs[int(self.v[1])-1].pos
        c = vs[int(self.v[2])-1].pos
        AB = np.array( [b[0]-a[0], b[1]-a[1], b[2]-a[2]] )
        AC = np.array( [c[0]-a[0], c[1]-a[1], c[2]-a[2]] )
        normal = np.cross(AB, AC)
        return normal / np.linalg.norm(normal)

class vFacesDict:
    def __init__(self, size, fs):
        logging.debug('创建每个顶点对应的面片集合...')
        logging.debug(f'size: {size}')
        self.vFaces = {}
        # regard str of the index as a key
        for i in range(1, size+1):
            self.vFaces[str(i)] = set()
        for f in fs:
            assert isinstance(f, Face)
            self.vFaces[str(f.v[0])].add(f)
            self.vFaces[str(f.v[1])].add(f)
            self.vFaces[str(f.v[2])].add(f)
        logging.debug(f'len of vFacesDict: {len(self.vFaces)}')
        logging.debug('顶点对应的面片集合创建成功...')

def Load_OBJ(filePath):
    vs = [] # vertex
    fs = [] # face's vertex
    vns = [] # n of vertex
    logging.debug('-----------------------------------------')
    logging.debug(f'正在加载模型文件: {filePath} ...')

    with open(filePath, 'r') as file:
        logging.debug('原模型打开成功...')
        for line in file.readlines():
            if line[0] == '#' or line[0] == '\n': continue
            line = line.rstrip('\n')
            datas = line.split(' ')
            # vertex
            if 'v' == datas[0]:
                v = Point(np.array( [float(datas[1]), float(datas[2]), float(datas[3])] ) ) # no need to save order
                vs.append(v)
                continue
            # n of vertex
            elif 'vn' == datas[0]:
                vns.append( np.array( [float(datas[1]), float(datas[2]), float(datas[3])] ) )
                continue
            # vt, ignore right now
            elif 'vt' == datas[0]: continue
            # f
            elif 'f' == datas[0]:
                buffer = [[], [], []] # 0->v; 1->vt; 2->vn
                for i in range(1, len(datas)):
                    temp = datas[i].split('/')
                    for ii in range(len(temp)):
                        buffer[ii].append(temp[ii])
                fs.append(Face( buffer[0], buffer[1], buffer[2] ))
            else: continue
        logging.debug('原模型加载完毕...')
        return vs, fs, vns
    
class MeshSubdivision:
    def __init__(self, vs=[], vns=[], fs=[], filepath=None):
        logging.debug('---------------------------------------------------')
        logging.debug('MeshSubdivision init...')
        self.sourcePath = filepath # set data path, for saving
        self.vs = vs
        self.vns = vns
        self.fs = fs
        self.vfs = vFacesDict(len(self.vs), self.fs)
        # list of new point
        self.newPoints = [] 
        # find the new point by a line which composes by 2 vertex
        self.newPointLineDict = {} 
        # list of new Faces
        self.newFaces = []
        logging.debug('MeshSubdivision is inited...')
    
    def getAdjacencies(self, ver:int)->set:
        """
        找到所有与点ver直接相连的点的集合
        ver 是点在self.vs的索引, 应为int
        """
        verKey = str(ver)
        vfaces = list(self.vfs.vFaces[verKey])
        adjacencies = set()
        # 点ver所在的每个三角 
        for f in vfaces:
            # 三角f的每个点
            for v in f.v:
                adjacencies.add(v)
        return adjacencies - set([ver])
    
    def commonFacesOf2Ver(self, v1:int, v2:int)->set:
        """ 求有两个相同顶点(v1, v2)的两个面 """
        t1 = self.vfs.vFaces[str(v1)]
        t2 = self.vfs.vFaces[str(v2)]
        res = set(t1) & set(t2) # union
        return res
    
    def updateVns(self)->None:
        """ 更新顶点法向量"""
        if len(self.vns)==0: return # 原模型无顶点法向量，无需update
        logging.debug('-----------------------------------------')
        logging.debug('正在更新顶点法向量')
        # 前置条件检查
        assert len(self.newFaces) != 0 # 保证新的面已经计算
        assert len(self.vs) > len(self.newPoints) #新旧顶点已经合并
        # update
        vns = []
        totalVerSize = len(self.vs) #新旧顶点总数
        self.vfs = vFacesDict(totalVerSize, self.newFaces)

        for i in range(1, totalVerSize+1):
            facesOfVer = list(self.vfs.vFaces[str(i)])
            # 孤立点
            if len(facesOfVer) == 0:
                vns.append(np.array( [0, 0, 0] ))
                continue
            Normals = [f.getNormal(self.vs) for f in facesOfVer ]
            res = sum(Normals) / len(facesOfVer)
            vns.append(res)
        self.vns = vns
        logging.debug('顶点法向量更新完毕...')

    def getNewPoints(self)->None:
        """ 计算新增点 """
        count = len(self.vs) + 1 #新增点的起始排序
        logging.debug('---------------------------------------')
        logging.debug('正在计算新增点...')
        for i in range(1, len(self.vs)+1):
            for verInAdj in self.getAdjacencies(i):
                # logging.debug(f'i={i}, ver={verInAdj}')
                # assert isinstance(verInAdj, int) ################
                t = [i, verInAdj]
                t.sort()
                # logging.debug(f't: {t}')
                t = str(t[0])+str(t[1]) #to str, as a key of dict
                if self.newPointLineDict.get(t, 'KNF') != 'KNF': continue # found, so it has been updated
                self.newPointLineDict[t] = Point(self.getPosBy2Ver(i, verInAdj), count)
                count+=1
        logging.debug(f'count: {count}')
        for item in self.newPointLineDict.keys():
            self.newPoints.append(self.newPointLineDict[item])
        logging.debug('新增点计算完成...')
        logging.debug(f'边数： {len(self.newPoints)}')

    def newPointBy2Ver(self, v1:int, v2:int)->Point:
        """ 根据两点, 返回该边的新增点"""
        t = [v1, v2]
        t.sort()
        t = str(t[0])+str(t[1])
        return self.newPointLineDict[t]

    def getPosBy2Ver(self, v1:int, v2:int)->np.ndarray:
        """ 根据两点, 计算相应的新增点坐标 """
        commonFaces = list(self.commonFacesOf2Ver(v1, v2))
        size = len(commonFaces)
        assert  size == 2 or size == 1
        # len() == 1, 该边为边界
        A = self.vs[v1 - 1].pos
        B = self.vs[v2 - 1].pos
        if size == 1:
            
            points = commonFaces[0].v
            farPoint = list(set(points) - set([v1, v2]))
            assert len(farPoint) == 1
            C = self.vs[farPoint[0] - 1].pos
            res = 7/15*(A+B) + 1/15*C
        # len() == 2, 非边界
        else:
            farPoint = list( set(commonFaces[0].v) & set(commonFaces[1].v))
            assert len(farPoint) == 2
            C = self.vs[farPoint[0] - 1].pos
            D = self.vs[farPoint[1] - 1].pos
            res = 3/8*(A + B) + 1/8*(C + D)
        return np.array(res)

    def updateOldPoint(self):
        """ 更新旧顶点并合并新旧顶点 """
        logging.debug('------------------------------------')
        logging.debug('___ in function updateOldPoint ___')
        logging.debug('正在合并新旧顶点...')
        oldLength = len(self.vs) #旧顶点数
        self.newPoints.sort()
        for item in self.newPoints:
            self.vs.append(item)
        logging.debug('新旧顶点合并完成...')
        
        logging.debug('正在更新旧顶点...')
        # oldV = []
        for i in range(1, oldLength+1):
            origin = self.vs[i-1].pos
            u= 0.0
            adj = list(self.getAdjacencies(i))
            adj = [self.newPointBy2Ver(i, ii).pos for ii in adj]
            n = len(adj)
            if n==0: continue
            elif n == 3:u = 3/16
            else: u = 1/(2*n)
            neighbor_position_sum = np.array([0.0, 0.0, 0.0])
            for item in adj:
                neighbor_position_sum = neighbor_position_sum + item
            self.vs[i-1].pos = ( (1-n*u)*origin + u*neighbor_position_sum)
        logging.debug('旧顶点更新完成...')
        logging.debug('顶点操作完成...')

    def getNewFaces(self):
        """ 计算新增面片 """
        logging.debug('----------------------------------')
        logging.debug('正在计算新增面...')
        for f in self.fs:
            temp = set()
            Vset = set(f.v) # 面f的顶点集合
            newPoints = set()
            for ver in f.v:
                other = list(Vset - set([ver]))
                assert len(other) == 2

                n1 = self.newPointBy2Ver(ver, other[0]).order
                n2 = self.newPointBy2Ver(ver, other[1]).order
                newPoints.add(n1)
                newPoints.add(n2)
                temp.add(Face([ver, n1, n2]))
            newPoints = list(newPoints)
            self.newFaces.append( Face([newPoints[0], newPoints[1], newPoints[2]]))
            for item in list(temp):
                self.newFaces.append(item)
        logging.debug('新增面计算完成...')

    def showDatas(self):
        logging.debug('--------------showing datas-------------------------')
        logging.debug(f'len of vs: { len(self.vs)}')
        logging.debug(f'len of vns: {len(self.vns)}')
        logging.debug(f'len of fs: {len(self.fs)}')
        logging.debug(f'len of nps: {len(self.newPoints)}')
        logging.debug(f'len of nfs: {len(self.newFaces)}')

    def storeAsOBJ(self, times):
        assert isinstance(times, int)
        assert self.sourcePath is not None
        targetPath = './results/' + f'sub_{times}_' + self.sourcePath
        logging.debug('---------------------------------------------------')
        logging.debug('正在存储为.obj文件...')
        logging.debug(f'目标路径为: {targetPath}')

        # 保证新旧顶点已经合并
        assert len(self.vs) > len(self.newPoints)
        # 新顶点应经过排序
        saveVns = len(self.vns) != 0

        with open(targetPath, 'w') as res:
            res.write('# here is result.obj\n')
            res.write(f'# verteice: {len(self.vs)}\n')
            res.write(f'# faces : {len(self.newFaces)}\n\n\n')
            # save points
            logging.debug('正在存储旧、新顶点...')
            for i in range(len(self.vs)):
                res.write(str(self.vs[i]) + '\n')
                if not saveVns: continue
                item = self.vns[i]
                temp = ['vn', f'{item[0]:.6f}', f'{item[1]:.6f}', f'{item[2]:.6f}']
                temp = ' '.join(temp)
                res.write(temp+'\n')
            # save faces
            logging.debug('正在存储面...')
            for item in self.newFaces:
                if not saveVns: item.vn = None
                res.write(str(item)+'\n')
        logging.debug('.obj文件已写入...')

    def runToObj(self, times):
        assert isinstance(times, int)
        self.getNewPoints()
        self.getNewFaces()
        self.updateOldPoint()
        self.updateVns()
        self.storeAsOBJ(times)

    def runTimes(self):
        self.getNewPoints()
        self.getNewFaces()
        self.updateOldPoint()

        # update
        self.fs, self.newFaces, self.newPoints = self.newFaces, [], []
        self.vfs = vFacesDict(len(self.vs), self.fs)
        self.newPointLineDict.clear()
    
def Subdivision(rootPath, filePath, times = 1):
    logging.debug('########################################')
    vs, fs, vns = Load_OBJ(rootPath + filePath)
    m = MeshSubdivision(vs, vns, fs, filePath)
    del vs, fs, vns
    m.showDatas()
    for i in range(times-1):
        logging.debug('########################################')
        logging.debug(f'Now times: {i}..........')
        m.runTimes()
        m.showDatas()
    logging.debug('___________________________________________________________')
    logging.debug('__________________________Final____________________________')
    m.runToObj(times)
    m.showDatas()

if __name__ == '__main__':
    rootPath = './HW2_datasets/Mesh/'
    filePath = 'Horse.obj'
    # rootPath = './Test/'
    # filePath = 'temp.obj'
    Subdivision(rootPath, filePath)