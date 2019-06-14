#include "pch.h"
#include "predefs.h"

// 开启时会进行一些运行时 assert 检查
constexpr bool DEBUG_MODE = true;

// 万有引力常数
constexpr double GRAVITY_CONST = 6.672e-11;

constexpr double THETA = 0.5, MASS_RATIO = 1000;

// 获取 pos 在 (begin, end) 的区域内对应的
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
    double gravityScalar = GRAVITY_CONST * A.mass * B.mass / distSquare;

    Vector2 gravity = v1to2.unit() * gravityScalar;
    return gravity;
}

class BHPool {
public:
    BHPool(int16_t bodyNum)
        : nodeIndexOfBody(std::make_unique<int16_t[]>(bodyNum)) {
        // 建立根节点
        nodes.push_back(BHNode {});

        // 根节点尽管没有子树，但特殊对待，看作内部节点
        nodes[0].internal = true;
    }

    void insert(const Body& body) {
        insert_on_node(body, 0, Vector2 { 0, 0 }, Vector2 { 1, 1 });
    }

    BHNode& operator[](const size_t nodeId) {
        return nodes[nodeId];
    }

    // 返回 Body 受到其它天体的引力之合力。
    Vector2 gravity(const Body& body) {
		return this->gravity_from_node(body, 0, 1.0);
    }

private:
    std::vector<BHNode> nodes;
    std::unique_ptr<int16_t[]> nodeIndexOfBody;

    int16_t add() {
        nodes.push_back(BHNode {});
        return static_cast<int16_t>(nodes.size() - 1);
    }

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
            Body bodyInNote = nodes[node].vaBody;

            // 插入（node, bodyInNode）
            insert_on_node(bodyInNote, node, begin, end);

            // 插入（node，body）
            insert_on_node(body, node, begin, end);
        }

        // 更新当前 node 的 vBody
        update_mass_vpos(node);
    }

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

    Vector2 gravity_from_node(const Body& body, const int16_t node, const double areaLength) {
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

    bool treat_as_single_body(const Body& body, int16_t node, double areaLength) {
        double dist = (nodes[node].vaBody.pos - body.pos).len();
        return (areaLength / dist) < THETA;
    }
};

void gen_body(Body& vBody) {
    std::mt19937_64 rnd { static_cast<uint64_t>(std::time(nullptr)) };
    std::uniform_real_distribution<> disPos(0.0, 1.0);

    vBody.pos.x = disPos(rnd);
    vBody.pos.y = disPos(rnd);

    std::uniform_real_distribution<> disMass(1.0, MASS_RATIO);
    vBody.mass = disMass(rnd);
}

void worker(const Settings& setting, int rank, int size) {
    // 1. 生成天体
    std::unique_ptr<Body[]> all_bodies;
    if (rank == 0) {
        all_bodies.reset(new Body[setting.bodyNum]);

        for (int i = 0; i < setting.bodyNum; i++) {
            gen_body(all_bodies[i]);
        }
    }

    // 2. 广播天体信息
    size_t bodyNumEach = setting.bodyNum / size;
    auto myBodies = std::make_unique<Body[]>(bodyNumEach);
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            size_t itsBody = i * bodyNumEach;
            MPI_Send(&myBodies[itsBody],
                static_cast<int>(bodyNumEach * sizeof(Body)), MPI_BYTE,
                i, 0, MPI_COMM_WORLD);
        }
        std::copy(all_bodies.get(), all_bodies.get() + bodyNumEach, myBodies.get());
    } else {
        MPI_Recv(&myBodies, static_cast<int>(bodyNumEach * sizeof(Body)), MPI_BYTE,
            0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }



    //// 3. 建树
    //BHPool* bhPool = nullptr;

    //if (rank == 0) {
    //    BHPool* bhPool = new BHPool(setting.bodyNum);
    //    for (int i = 0; i < setting.bodyNum; i++) {
    //        bhPool->insert(all_bodies[i]);
    //    }
    //}
}

int mpi_main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Settings setting = { 0 };
    if (rank == 0) {
        std::cout << "天体数量：";
        std::cin >> setting.bodyNum;

        // 要保留一定空间
        if (setting.bodyNum >= std::numeric_limits<int16_t>::max()) {
            std::cout << "ERROR：天体数量太多。\n退出……\n"
                      << std::endl;
            return -1;
        }

        // 天体必须能平均分给每个节点
        if (setting.bodyNum % size != 0) {
            std::cout << "ERROR：天体不能被 size 整除。\n退出……\n"
                      << std::endl;
            return -1;
        }
    }
    MPI_Bcast(&setting, sizeof(setting), MPI_BYTE, 0, MPI_COMM_WORLD);

    std::cout << "[" << rank << "] bodyNum = " << setting.bodyNum;
    worker(setting, rank, size);

    MPI_Finalize();
    return 0;
}

void bhpool_test() {
    BHPool& pool = *new BHPool(5);
    std::vector<Body> bodies;

    bodies.push_back({ { 0.25, 0.25 }, 1 });
    bodies.push_back({ { 0.9, 0.37 }, 1 });
    bodies.push_back({ { 0.8, 0.37 }, 1 });
    bodies.push_back({ { 0.8, 0.12 }, 1 });
    bodies.push_back({ { 0.75, 0.75 }, 1 });

    for (Body& body : bodies) {
        pool.insert(body);
    }

    Vector2 totalPos;
    double totalMass = 0.0;
    for (Body& body : bodies) {
        totalPos += body.pos * body.mass;
        totalMass += body.mass;
    }
    totalPos /= totalMass;
    std::cout << totalPos.x << ", " << totalPos.y << " | " << totalMass << std::endl;
    std::cout << pool[0].vaBody.pos.x << ", " << pool[0].vaBody.pos.y << " | " << pool[0].vaBody.mass;
}

int main(int argc, char* argv[]) {
    // return mpi_main(argc, argv);
}
