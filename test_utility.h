#include <string>
#include <vector>
#include "./nlohmann/json.hpp"

using json = nlohmann::json;

class TestDescription{
public:
    TestDescription(std::string test_name, bool vertex_labelled, bool directed, bool connected, bool edge_labelled, float timeout, std::string node_heuristic):
    test_name(test_name), vertex_labelled(vertex_labelled), directed(directed), connected(connected), 
    edge_labelled(edge_labelled), timeout(timeout), node_heuristic(node_heuristic){};

    TestDescription(){}
    std::string test_name;
    bool vertex_labelled;
    bool directed;
    bool connected;
    bool edge_labelled;
    float timeout;
    std::string node_heuristic;
};

class Test{
public:
    Test(TestDescription td):td(td), total_time(0), recursions(0){}
    Test(){}
    void to_string();
    TestDescription td;
    float total_time;
    unsigned long long recursions;
    std::vector<std::tuple<int, unsigned long long, double>> milestones;
    std::vector<std::pair<int, int>> solution;
};

inline void to_json(json& j, const TestDescription& t){
    j = json{{"test_name", t.test_name}, {"vertex_labelled", t.vertex_labelled}, {"directed", t.directed}, {"connected", t.connected},
        {"edge_labelled", t.edge_labelled}, {"timeout", t.timeout}, {"node_heuristic", t.node_heuristic}};
}

inline void to_json(json& j, const Test& t){
    j = json{{"milestones", t.milestones}, {"solution", t.solution}, {"test_description", t.td}, {"recursions", t.recursions},
        {"total_time", t.total_time}};
}

void save_json(Test& t, std::string save_folder);
