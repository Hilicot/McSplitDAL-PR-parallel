#include "test_utility.h"
#include <iostream>

#include <fstream>

void Test::to_string(){
    std::cout << std::endl << "TEST DESCRIPTION" << std::endl;
    std::cout << "Test Name: " << this->td.test_name << " vertex_labelled: " << this->td.vertex_labelled 
        << " directed: " << this->td.directed << " connected: " << this->td.connected << " edge_labelled: " << this->td.edge_labelled
        << " timeout: " << this->td.timeout << " node_heuristic: " << this->td.node_heuristic << std::endl << std::endl;

    std::cout << "Total Number of Recursions: " << this->recursions << std::endl;
    std::cout << "Total Elapsed Time: " << this->total_time << std::endl;
    std::cout << "MILESTONES" << std::endl;
    for(auto& m : this->milestones){
        std::cout << "Incumbent Size: " << std::get<0>(m) << " Recursions: " << std::get<1>(m) << " Time: " << std::get<2>(m) << std::endl;
    }

    std::cout << "SOLUTION" << std::endl;
    std::cout << "Solution size: " << this->solution.size() << std::endl;
    for(auto& m : this->solution){
        std::cout << "(" << m.first << " -> " << m.second << ") ";
    }
}

void save_json(Test& t, std::string save_folder){
    json jsonfile = t;

    std::ofstream file(save_folder + "/" + t.td.test_name  + ".json");
    file << jsonfile;
 }