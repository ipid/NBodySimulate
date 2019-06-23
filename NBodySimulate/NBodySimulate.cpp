#include "pch.h"
#include "predefs.h"

#ifdef _DEBUG
constexpr bool DEBUG_MODE = true;
#else
constexpr bool DEBUG_MODE = false;
#endif

// 打开时会运行测试代码
constexpr bool TEST_MODE = false;

// 打开此开关可以只测试性能，不输出天体结果
constexpr bool BENCHMARK_ONLY_AND_DO_NOT_PRINT_BODY_RESULT = false;

constexpr double THETA = 0.25, MASS_RATIO = 1e4;
constexpr double MOTION_DELTA_TIME = 0.1, INITIAL_VELOCITY = 0.02;
constexpr int INDENT = 2;

enum Messages : int {
    InitialBodyInfo = 0,
    NodeNumOfCurrentIter,
    RawNodesOfCurrentIter,
    SlaveBodyStatus
};

template <class T>
constexpr int int_sizeof(const T& elem, size_t num) {
    return static_cast<int>(num * sizeof(T));
}

template <class T>
constexpr int int_sizeof(size_t num) {
    return static_cast<int>(num * sizeof(T));
}

// 计算传入的 body 是否在大边界内
// 超出边界范围的 body 不应参与计算
bool should_body_be_calculated(const Body& body) {
    return (body.pos.x > .0 && body.pos.y > .0 && body.pos.x < 1.0 && body.pos.y < 1.0);
}

// 打印 body 数组中所有天体的位置信息
// 超出边界的 body 位置用“-”表示
bool print_body_status(const Body body[], const size_t bodyNum) {
	bool isAnyBodyPrinted = false;

    for (size_t i = 0; i < bodyNum; i++) {
        if (should_body_be_calculated(body[i])) {
            isAnyBodyPrinted = true;
            std::cout << body[i].pos.x << ' ' << body[i].pos.y << ' ';
        } else {
            std::cout << "- - ";
        }
    }
    std::cout << std::endl;
    return isAnyBodyPrinted;
}

void put_spaces(size_t spaceNum) {
    for (size_t i = 0; i < spaceNum; i++) {
        std::cout.put(' ');
    }
}

// 获取 pos 在 (begin, end) 的区域内对应的象限
Quadrant which_quadrant(const Vector2& pos,
    const Vector2& begin, const Vector2& end) {

    Vector2 middle = (begin + end) / 2;

    if constexpr (DEBUG_MODE) {
        if (!end.at_rd_of(begin)) {
            throw "你妈的，为什么";
        }
    }

    if (pos.x < middle.x && pos.y < middle.y) {
        return Quadrant::LU;
    } else if (pos.x < middle.x && pos.y >= middle.y) {
        return Quadrant::LD;
    } else if (pos.x >= middle.x && pos.y < middle.y) {
        return Quadrant::RU;
    } else if (pos.x >= middle.x && pos.y >= middle.y) {
        return Quadrant::RD;
    }

    if constexpr (DEBUG_MODE) {
        throw "你妈的，为什么";
    }
}

// 输入区域和象限 quadr，该函数将区域分为四个子区域，返回对应 quadr 象限的子区域
std::pair<Vector2, Vector2> border_of_quadrant(const Vector2& begin,
    const Vector2& end, Quadrant quadr) {

    if constexpr (DEBUG_MODE) {
        if (!end.at_rd_of(begin)) {
            throw "你妈的，为什么";
        }
    }

    Vector2 middle = (begin + end) / 2;

    switch (quadr) {
    case Quadrant::RU:
        return { { middle.x, begin.y }, { end.x, middle.y } };
    case Quadrant::RD:
        return { middle, end };
    case Quadrant::LU:
        return { begin, middle };
    case Quadrant::LD:
        return { { begin.x, middle.y }, { middle.x, end.y } };
    }
}

// 返回 A、B 之间的引力，方向为从 A 到 B
Vector2 newton_gravity_AtoB(const Body& A, const Body& B) {
    Vector2 v1to2 = B.pos - A.pos;

    // GRAVITY_CONST * body1.mass * body2.mass / (body2 - body1) ^ 2
    double distSquare = v1to2.len_square();

    // 如果两个点正好重合
    if (distSquare <= 0.0) {
        // 就让这两个点合为一体，或者通过惯性分开
        return { 0, 0 };
    }

    double gravityScalar = GRAVITY_CONST * A.mass * B.mass / distSquare;

    Vector2 gravity = v1to2.unit() * gravityScalar;
    return gravity;
}

// 使用已建立的Barnes-Hut树计算天体引力
class BHPoolUtility {
public:
    BHPoolUtility(BHNode nodes[], const size_t nodeNum)
        : nodes(nodes)
        , nodeNum(nodeNum) {
    }

    BHPoolUtility(const BHPoolUtility& other)
        : nodes(other.nodes)
        , nodeNum(other.nodeNum) {
    }

    // 返回 Body 受到其它天体的引力之合力。
    Vector2 gravity(const Body& body) const {
        if (!should_body_be_calculated(body)) {
            return { .0, .0 };
        }
        return this->gravity_from_node(body, 0, 1.0);
    }

    void display_tree() const {
        this->display_tree_from_node(0, 0);
        std::cout.flush();
    }

private:
    BHNode* nodes;
    size_t nodeNum;

	// 计算节点 node 对 body 产生的引力
    Vector2 gravity_from_node(const Body& body, const int16_t node, const double areaLength) const {
        if (nodes[node].vaBody == body) {
            return { .0, .0 };
        }

        // 如果为外部节点或者 treat_as_single_body
        if ((!nodes[node].internal) || treat_as_single_body(body, node, areaLength)) {
            return newton_gravity_AtoB(body, nodes[node].vaBody);
        }

        Vector2 total;
        for (int i = 0; i < 4; i++) {
            if (nodes[node].child[i] < 0) {
                continue;
            }

            total += gravity_from_node(body, nodes[node].child[i], areaLength / 2);
        }

        return total;
    }

	// node 在与 body 相互作用时，是否应被视为一整个物体
    bool treat_as_single_body(const Body& body, int16_t node, double areaLength) const {
        double dist = (nodes[node].vaBody.pos - body.pos).len();
        return (areaLength / dist) < THETA;
    }

	// DEBUG: 打印当前 node 树
    void display_tree_from_node(int16_t node, int indent) const {
        put_spaces(indent);
        if (node < 0) {
            std::cout << "- ,\n";
            return;
        }

        std::cout << "( " << nodes[node].vaBody << ", {\n";
        for (size_t i = 0; i < 4; i++) {
            display_tree_from_node(nodes[node].child[i], indent + INDENT);
        }

        put_spaces(indent);
        std::cout << "})\n";
    }
};

// 建立Barnes-Hut树
class BHPool {
public:
    BHPool() {
        // 建立根节点
        nodes.push_back(BHNode {});

        // 根节点尽管没有子树，但特殊对待，看作内部节点
        nodes[0].internal = true;
    }

	// 往树中插入天体信息
    void insert(const Body& body) {
        if (should_body_be_calculated(body)) {
            insert_on_node(body, 0, Vector2 { 0, 0 }, Vector2 { 1, 1 });
        }
    }

	// 直接访问 nodeId 对应的 BH 树节点
    BHNode& operator[](const size_t nodeId) {
        return nodes[nodeId];
    }

	// 复制一份内部的节点并返回
    std::unique_ptr<BHNode[]> export_raw_nodes_copy() {
        auto rawNodes = std::make_unique<BHNode[]>(nodes.size());
        std::copy(nodes.begin(), nodes.end(), rawNodes.get());
        return rawNodes;
    }

	// 返回 BH 树节点数量
    size_t node_num() {
        return nodes.size();
    }

private:
    std::vector<BHNode> nodes;

	// 动态创建 BH 树节点
    int16_t add() {
        nodes.push_back(BHNode {});
        return static_cast<int16_t>(nodes.size() - 1);
    }

	// 在 node 上插入 body
    void insert_on_node(const Body& body, const int16_t node,
        const Vector2& begin, const Vector2& end) {

        // quadrOfNode = body 在 bound 之内所在象限
        Quadrant quadrBodyIn = which_quadrant(body.pos, begin, end);
        size_t quadrInt = static_cast<size_t>(quadrBodyIn);

        // 计算 quadrBodyIn 的边界
        auto [newBegin, newEnd] = border_of_quadrant(begin, end, quadrBodyIn);

        // 如果当前 node 为内部节点
        if (nodes[node].internal) {

            // 如果 child[quadrOfNode] 上没有节点
            if (nodes[node].child[quadrInt] < 0) {

                // 建立新节点，然后赋值给 child[quadrOfNode]
                int16_t newNode = this->add();
                nodes[node].child[quadrInt] = newNode;

                // 将 body 放入新节点
                nodes[newNode].vaBody = body;
            } else {

                // 插入（child[quadrOfNode]，body）
                this->insert_on_node(body, nodes[node].child[quadrInt], newBegin, newEnd);
            }

        } else {
            // 将 node 转化（即标记）为内部节点
            nodes[node].internal = true;

            // 取出 node 内的 body -> bodyInNode（此时 vaBody 为实际 body）
            Body bodyInNode = nodes[node].vaBody;

            // 插入（node, bodyInNode）
            insert_on_node(bodyInNode, node, begin, end);

            // 插入（node，body）
            insert_on_node(body, node, begin, end);
        }

        // 更新当前 node 的 vBody
        update_mass_vpos(node);
    }

	// 更新当前 node 上的虚拟天体信息
    void update_mass_vpos(const int16_t node) {
        Vector2 totalPosition;
        double totalMass = 0;

        for (int i = 0; i < 4; i++) {
            int16_t& child = nodes[node].child[i];

            if (child >= 0) {
                totalPosition += nodes[child].vaBody.pos * nodes[child].vaBody.mass;
                totalMass += nodes[child].vaBody.mass;
            }
        }

        nodes[node].vaBody.pos = totalPosition / totalMass;
        nodes[node].vaBody.mass = totalMass;
    }
};

// 随机生成天体
// 随机生成位置、质量、初始速度
void generate_bodies(Body bodies[], size_t bodyNum) {
    std::mt19937_64 rnd { static_cast<uint64_t>(std::time(nullptr)) };

    // 位置、质量、初始速度角度
    std::uniform_real_distribution<> disPos(0.0, 1.0),
        disMass(1.0, MASS_RATIO), disAngle(0, PI_CONST * 2);

    for (int i = 0; i < bodyNum; i++) {
        Body& vBody = bodies[i];

        vBody.pos.x = disPos(rnd);
        vBody.pos.y = disPos(rnd);
        vBody.mass = disMass(rnd);

        // 初始速度大小固定，方向随机产生（单位为弧度）
        double velocityAngle = disAngle(rnd);
        vBody.v = { INITIAL_VELOCITY * std::cos(velocityAngle),
            INITIAL_VELOCITY * std::sin(velocityAngle) };
    }
}

// 使用 BHUtility 来更新 myBodies 数组里的 Body 状态（有副作用）
// delta_time 单位为“秒”
void update_body_status(const BHPoolUtility& util,
    Body myBodies[], size_t bodyNumEach, double delta_time = MOTION_DELTA_TIME) {

    for (int i = 0; i < bodyNumEach; i++) {
        if (!should_body_be_calculated(myBodies[i])) {
            continue;
        }

        Vector2 gravity = util.gravity(myBodies[i]);
        myBodies[i].move_under_force(gravity, delta_time);
    }
}

/*
	[root]：表示 rank 为 0 的进程要做的事
	[slave]：表示 rank 不为 0 的进程都要做的事
	[worker]：[root] + [slave] 全体进程都要做的事
*/
void worker(const Settings& setting, int rank, int size) {
    if (rank == 0) {
        std::cout << "\n[BEGIN]" << std::endl;	
	}
    std::unique_ptr<Body[]> allBodies;

    // 1. [root] 生成天体，并打印初始天体信息
    if (rank == 0) {
        allBodies = std::make_unique<Body[]>(setting.bodyNum);
        generate_bodies(allBodies.get(), setting.bodyNum);

        if constexpr (!BENCHMARK_ONLY_AND_DO_NOT_PRINT_BODY_RESULT) {
            print_body_status(allBodies.get(), setting.bodyNum);
        } else {
            std::cout << "Begin: " << MPI_Wtime() << std::endl;
        }
    }

    // 2. [root] 逐个发送天体信息 / [slave] 接收天体信息
    const size_t bodyNumEach = setting.bodyNum / size;
    auto myBodies = std::make_unique<Body[]>(bodyNumEach);
    if (rank == 0) {
        for (int slave = 1; slave < size; slave++) {
            // 将属于第 slave 号 worker 的天体发给他

            size_t firstBodyOfSlave = slave * bodyNumEach;
            MPI_Send(&allBodies[firstBodyOfSlave],
                int_sizeof<Body>(bodyNumEach), MPI_BYTE,
                slave, InitialBodyInfo, MPI_COMM_WORLD);
        }
        // 将属于根节点的节点复制到 myBodies 中
        std::copy(allBodies.get(), allBodies.get() + bodyNumEach, myBodies.get());
    } else {
        // 接收属于当前节点的天体信息
        MPI_Recv(myBodies.get(), int_sizeof<Body>(bodyNumEach), MPI_BYTE,
            0, InitialBodyInfo, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // 循环进行 3 - 6 的操作
    for (size_t currentIter = 0; currentIter < setting.iterations; currentIter++) {

        // 3. [root] 建立所有天体的 Barnes-Hut 树
        std::unique_ptr<BHNode[]> rawNodes;
        size_t nodeNum = std::numeric_limits<size_t>::max();

        if (rank == 0) {
            BHPool pool;
            for (int i = 0; i < setting.bodyNum; i++) {
                pool.insert(allBodies[i]);
            }

            nodeNum = pool.node_num();
            MPI_Bcast(&nodeNum, sizeof(nodeNum), MPI_BYTE, 0, MPI_COMM_WORLD);

            rawNodes = pool.export_raw_nodes_copy();
            MPI_Bcast(rawNodes.get(), int_sizeof<BHNode>(nodeNum), MPI_BYTE, 0, MPI_COMM_WORLD);
        } else {
            MPI_Bcast(&nodeNum, sizeof(nodeNum), MPI_BYTE, 0, MPI_COMM_WORLD);

            rawNodes = std::make_unique<BHNode[]>(nodeNum);
            MPI_Bcast(rawNodes.get(), int_sizeof<BHNode>(nodeNum), MPI_BYTE, 0, MPI_COMM_WORLD);
        }

        BHPoolUtility utility(rawNodes.get(), nodeNum);

        // 4. [worker] 计算每个天体的下一 DELTA_TIME 状态
        update_body_status(utility, myBodies.get(), bodyNumEach);

        // 5. [root] 将天体状态存入 allBodies、接收天体状态 / [slave] 将天体状态发给 root
        if (rank == 0) {
            std::copy(&myBodies[0], myBodies.get() + bodyNumEach, &allBodies[0]);

            for (int slave = 1; slave < size; slave++) {
                size_t firstBodyOfSlave = slave * bodyNumEach;

                MPI_Recv(&allBodies[firstBodyOfSlave], int_sizeof<Body>(bodyNumEach), MPI_BYTE,
                    slave, SlaveBodyStatus, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            MPI_Send(&myBodies[0], int_sizeof<Body>(bodyNumEach), MPI_BYTE, 0, SlaveBodyStatus, MPI_COMM_WORLD);
        }

        // 6. [root] 输出天体当前状态
        if constexpr (!BENCHMARK_ONLY_AND_DO_NOT_PRINT_BODY_RESULT) {
            if (rank == 0) {
                print_body_status(allBodies.get(), setting.bodyNum);
            }
        }
    }
    if (rank == 0) {
        std::cout << "[END]" << std::endl;
    }

    if constexpr (BENCHMARK_ONLY_AND_DO_NOT_PRINT_BODY_RESULT) {
        if (rank == 0) {
            std::cout << "Finish: " << MPI_Wtime() << std::endl;
        }
    }
}

int mpi_main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Settings setting;
    if (rank == 0) {
        // 1. 设置天体个数
        std::cout << "Body Num: ";
        std::cin >> setting.bodyNum;

        // 要保留一定空间
        if (setting.bodyNum >= std::numeric_limits<int16_t>::max()) {
            std::cout << "ERROR: Body num must be smaller than "
                      << std::numeric_limits<int16_t>::max()
                      << ".\nExit...\n"
                      << std::endl;
            return -1;
        }

        // 天体必须能平均分给每个节点
        if (setting.bodyNum % size != 0) {
            std::cout << "ERROR: Body num can not be divided by size = "
                      << size << ".\nExit...\n"
                      << std::endl;
            return -1;
        }

        // 2. 设置迭代次数
        std::cout << "Iteration: ";
        std::cin >> setting.iterations;

        if (setting.iterations < 3) {
            std::cout << "ERROR: Iteration time must be more than 3 times.\nExit...\n"
                      << std::endl;
            return -1;
        }
    }
    // 广播基础设定
    MPI_Bcast(&setting, sizeof(setting), MPI_BYTE, 0, MPI_COMM_WORLD);

    // 初始化完毕，开始工作
    worker(setting, rank, size);

    MPI_Finalize();
    return 0;
}

// DEBUG: BHPool 的测试代码，测试树是否能正确建立
void bhpool_test() {
    BHPool pool;
    std::vector<Body> bodies;

    bodies.push_back({ { 0.25, 0.25 }, 1 });
    bodies.push_back({ { 0.9, 0.37 }, 1 });
    bodies.push_back({ { 0.8, 0.37 }, 1 });
    bodies.push_back({ { 0.8, 0.12 }, 1 });
    bodies.push_back({ { 0.75, 0.75 }, 1 });
    for (Body& b : bodies) {
        pool.insert(b);
    }

    Vector2 gravity0;
    for (int i = 1; i < 5; i++) {
        gravity0 += newton_gravity_AtoB(bodies[0], bodies[i]);
    }

    auto rawNodes = pool.export_raw_nodes_copy();
    BHPoolUtility utility { rawNodes.get(), pool.node_num() };

    Vector2 gravity0BH = utility.gravity(bodies[0]);
    std::cout << gravity0 << "\n"
              << gravity0BH << std::endl;

    utility.display_tree();
}

// DEBUG: 模拟一个双星系统，从而测试引力、BH 树、动量定理等代码
void test_doubleStarSystem_ResultCorrect() {
    std::vector<Body> bodies;
    bodies.push_back(Body { { 0.5, 0.446624 }, { -0.5, .0 }, 8e8 });
    bodies.push_back(Body { { 0.5, 0.553376 }, { 0.5, .0 }, 8e8 });

    for (int count = 0; count < 3000; count++) {

        BHPool pool;
        for (Body& b : bodies) {
            pool.insert(b);
        }

        auto rawNodes = pool.export_raw_nodes_copy();
        BHPoolUtility util { rawNodes.get(), 2 };

        std::cout << bodies[0].pos.x << " "
                  << bodies[0].pos.y << " "
                  << bodies[1].pos.x << " "
                  << bodies[1].pos.y << std::endl;
        update_body_status(util, bodies.data(), 2, 0.01);
    }
}

int main(int argc, char* argv[]) {
    // 尽可能高精度
    std::cout.precision(10);

    if constexpr (TEST_MODE) {
        // bhpool_test();
        test_doubleStarSystem_ResultCorrect();

        std::cout << "Memory leak: "
                  << (_CrtDumpMemoryLeaks() ? "Leaked" : "No leaking")
                  << std::endl;
        return 0;
    } else {
        return mpi_main(argc, argv);
    }
}
