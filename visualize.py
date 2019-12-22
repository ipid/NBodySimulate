import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# 读取模拟文件
def read_simulate_res():
    res = []

    with open('simulate-res.txt') as dataFile:
        content = dataFile.read()
        dataLines = content.split('\n')
        dataLinesIter = iter(dataLines)

        # 读到 [BEGIN] 位置
        while next(dataLinesIter) != '[BEGIN]':
            pass

        for data in dataLinesIter:
            if data == '[END]':
                break

            bodyLocs = data.strip().split(' ')

            resX, resY = [], []
            for i, elem in enumerate(bodyLocs):
                if elem == '-':
                    continue

                if i % 2 == 0:
                    resX.append(float(elem))
                else:
                    resY.append(float(elem))

            res.append((resX, resY))

    return res


# 初始化 matplotlib
fig, ax = plt.subplots()
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
line, = ax.plot([], [], '.b')

frameData = read_simulate_res()

def ani_func(i):
    line.set_data(frameData[i][0], frameData[i][1])
    return line,


def ani_init():
    return ani_func(0)


anim = FuncAnimation(fig, ani_func, init_func=ani_init, frames=len(frameData), interval=1)
plt.show()
