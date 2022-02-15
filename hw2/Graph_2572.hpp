#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
class Graph {
private:

    // HW0: YOUR CODE HERE
    // Use this space for declarations of important internal types you need
    // later in the Graph's definition.
    // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
    // code here. Just use the space if you need it.)

    // Declare struct for node's internal data
    struct internal_node;
    struct internal_edge;

public:

    //
    // PUBLIC TYPE DEFINITIONS
    //

    /** Type of node value. */
    using node_value_type = V;

    /** Type of edge value. */
    using edge_value_type = E;

    /** Type of this graph. */
    using graph_type = Graph<V,E>;

    /** Predeclaration of Node type. */
    class Node;
    /** Synonym for Node (following STL conventions). */
    using node_type = Node;

    /** Predeclaration of Edge type. */
    class Edge;
    /** Synonym for Edge (following STL conventions). */
    using edge_type = Edge;

    /** Type of node iterators, which iterate over all graph nodes. */
    class NodeIterator;
    /** Synonym for NodeIterator */
    using node_iterator = NodeIterator;

    /** Type of edge iterators, which iterate over all graph edges. */
    class EdgeIterator;
    /** Synonym for EdgeIterator */
    using edge_iterator = EdgeIterator;

    /** Type of incident iterators, which iterate incident edges to a node. */
    class IncidentIterator;
    /** Synonym for IncidentIterator */
    using incident_iterator = IncidentIterator;

    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {
        // HW0: YOUR CODE HERE
        node_counter = 0;
        edge_counter = 0;
    }

    /** Default destructor */
    ~Graph() = default;

    //
    // NODES
    //

    /** @class Graph::Node
     * @brief Class representing the graph's nodes.
     *
     * Node objects are used to access information about the Graph's nodes.
     */
    class Node : private totally_ordered<Node>{
    public:
        /** Construct an invalid node.
         *
         * Valid nodes are obtained from the Graph class, but it
         * is occasionally useful to declare an @i invalid node, and assign a
         * valid node to it later. For example:
         *
         * @code
         * Graph::node_type x;
         * if (...should pick the first node...)
         *   x = graph.node(0);
         * else
         *   x = some other node using a complicated calculation
         * do_something(x);
         * @endcode
         */
        Node() {
            // HW0: YOUR CODE HERE
            is_valid_ = false;
        }


        /** Return this node's modifiable position. */
        Point& position() {
            return fetch().pos;
        }

        /** Return this node's position. */
        const Point& position() const {
            // HW0: YOUR CODE HERE
            return fetch().pos;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            return fetch().idx_;
        }

        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS

        /** Return this node's value, of type node_value_type. */
        node_value_type& value() {
            return fetch().v_type;
        }

        /** Return this constant node's value, of type node_value_type. */
        const node_value_type& value() const{
            const internal_node& node = fetch();
            return node.v_type;
        }

        /** Return this node's number of neighbors. */
        size_type degree() const {
            return fetch().neighbors.size();
        }

        /** Return the beginning of node's incident iterator over neighbors. */
        incident_iterator edge_begin() const {
            return IncidentIterator(parent_graph_, index(),0,&fetch().neighbors);
        }

        /** Return the end of node's incident iterator over neighbors. */
        incident_iterator edge_end() const {
            return IncidentIterator(parent_graph_, index(),
                                    degree(),
                                    &fetch().neighbors);
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            // HW0: YOUR CODE HERE
            bool c1 = fetch().idx_ == n.index();
            bool c2 = parent_graph_ == n.parent_graph_;
            return (c1 && c2);
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node& n) const {
            // HW0: YOUR CODE HERE
            if (fetch().idx_ < n.index()) {
                return true;
            }
            return false;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects

        // Pointer to graph that node is in
        graph_type* parent_graph_;
        // Unique id for node
        size_type n_id_;
        // Flag for valid id
        bool is_valid_;

        /** Private Constructor for Graph, creates valid node. */
        Node(const graph_type* parent_graph, size_type n_id)
                : parent_graph_(const_cast<Graph*>(parent_graph)), n_id_(n_id){
            is_valid_= true;
        }

        /** Helper method to return current internal node object. */
        internal_node& fetch() const {
            // return the internal_node that has the same index as n_id_
            assert(parent_graph_->node_index_map_.count(n_id_) != 0);
            size_type node_index = parent_graph_->node_index_map_.at(n_id_);
            // check if node_index is a valid id
            if (node_index < parent_graph_->node_vec_.size()) {
                return parent_graph_->node_vec_[node_index];
            } else {
                assert(false);
            }
        }

    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        // HW0: YOUR CODE HERE
        //return node_vec_.size();
        // TODO: figure this out 
        return node_vec_.size();
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @param[in] val The new node's value type
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position,
                  const node_value_type& val = node_value_type()) {
        // HW0: YOUR CODE HERE
        // create a new unique id 
        size_type new_id = node_counter;
        node_counter++;

        // find what index it will have 
        size_type new_idx = node_vec_.size();

        // update internal node attributes
        internal_node add_internal;
        add_internal.pos = position;
        add_internal.n_id = new_id;
        add_internal.v_type = val;
        add_internal.neighbors = {};
        add_internal.idx_ = new_idx;

        // keep track of the unique id -> index mapping
        node_index_map_.insert({new_id, new_idx});
        node_vec_.push_back(add_internal); // add internal attributes to vector
        return Node(this, new_id);
    }

    /** Remove a node to the graph, returning 0 if the node isn't in the graph.
     * @param[in] n The node we want to remove
     * @post new num_nodes() == old num_nodes() - 1
     *
     * Complexity: O(num_nodes())
     */
    size_type remove_node(const Node& n) {
        // check if node is in the graph. if it isn't return 0
        if (n.index() >= num_nodes()){
            return 0;
        }  

        // keep track of the index, unique id, and internal data of the node we
        //     want to delete. 
        size_type idx = n.index();
        size_type id = node_vec_[idx].n_id;
        internal_node old_int = node_vec_[idx];

        // make a copy of the node's neighbors to loop through
        std::vector<size_type> nbs = old_int.neighbors;

        // for every neighbor, remove the edge between the node and neighbor
        for (unsigned i = 0; i < nbs.size(); i++){
            size_type n_idx = node_index_map_.at(nbs[i]);
            remove_edge(n, node(n_idx));
        }
        
        // swap the node vec to overwrite the node_internal we want to remove 
        node_vec_[idx] = node_vec_.back();
        node_vec_[idx].idx_ = idx;
        // update the node mapping
        node_index_map_.at(node_vec_[idx].n_id) = idx;

        // pop the node vec to remove repeated info
        node_vec_.pop_back();
        // remove the old mapping from old unique id
        node_index_map_.erase(id);
        return 1;
    }

    /** Remove a node to the graph, returning 0 if the node isn't in the graph.
     * @param[in] n The node we want to remove
     * @post new num_nodes() == old num_nodes() - 1
     *
     * Complexity: O(num_nodes())
     */
    node_iterator remove_node(node_iterator n_it) {
        // dereference iterator and call other remove node method
        size_type b = remove_node(*(n_it));
        if (b == 0){
            return n_it;
        } else {
            return node_begin();
        }
    }


    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        bool c1 = (n.index() < num_nodes());
        bool c2 = (n.parent_graph_ == this);
        return (c1 && c2);
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
        if (i < num_nodes()){
            size_type id = node_vec_[i].n_id;
            return Node(this,id);
        } else {
            assert(false);
        }
    }

    //
    // EDGES
    //

    /** @class Graph::Edge
     * @brief Class representing the graph's edges.
     *
     * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
     * are considered equal if they connect the same nodes, in either order.
     */
    class Edge : private totally_ordered<Edge>{
    public:
        /** Construct an invalid Edge. */
        Edge() {
            // HWO: YOUR CODE HERE
            is_valid_ = false;
        }

        /** size_type index() const {
            return fetch().idx_;
        }*/

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            size_type node1_id = fetch().an_edge_.first;
            size_type node1_idx = parent_graph_->node_index_map_.at(node1_id);
            return parent_graph_->node(node1_idx);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            size_type node2_id = fetch().an_edge_.second;
            size_type node2_idx = parent_graph_->node_index_map_.at(node2_id);
            return parent_graph_->node(node2_idx);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            // HW0: YOUR CODE HERE
            if (parent_graph_ != e.parent_graph_) {
                return false;
            }
            if ((node1() == e.node1()) && (node2() == e.node2())) {
                return true;
            } else if ((node2() == e.node1()) && (node1() == e.node2())) {
                return true;
            }
            return false;
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            // HW0: YOUR CODE HERE
            if (*this == e) {
                return false;
            } else if ((parent_graph_)!=e.parent_graph_) {
                return (parent_graph_ < e.parent_graph_);
            }else if ((node1() < e.node1()) && (node2() < e.node2())) {
                return true;
            } 
            return false;
        }

        /**
         * @brief Return the length of the edge. Length is calculated using the
         * L2-norm between the two node's positions. 
         * 
         * @return double >= 0
         */
        double length() const {
            Point a = node1().position();
            Point b = node2().position();
            return norm(b - a);
        }

        /**
         * Method to return the value of the edge
         */
        const edge_value_type& value() const {
            const internal_edge& edge = fetch();
            return edge.e_val;
        }

        /**
         * Non-const method to return the value of the edge. Note we can 
         * modify the edge's internal value with this method.
         */
        edge_value_type& value() {
            return fetch().e_val;
        }


    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members sand methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects

        graph_type* parent_graph_;
        size_type e_id_;
        bool is_valid_;


        /** Private constructor so Graph can create new edges */
        Edge(const graph_type* parent_graph, size_type e_id) :
                    parent_graph_(const_cast<Graph*>(parent_graph)), 
                    e_id_(e_id) {
                        is_valid_ = true;
                    }

        /** Helper method to update the current orientation of edge. 
         * @pre @a a and @a b must be between 0 and node_counter, 
         * i.e. they are unique node id's
        */
        void update_edge(size_type a, size_type b) {
            fetch().an_edge_.first = a;
            fetch().an_edge_.second = b;
        }

        /** Helper method to return current internal edge object. */
        internal_edge& fetch() const {
            // return the internal edge that has the same e_id as e_id_
            assert(parent_graph_->edge_index_map_.count(e_id_) != 0);
            size_type edge_index = parent_graph_->edge_index_map_.at(e_id_);
            if (edge_index < parent_graph_->edge_vec_.size()) {
                return parent_graph_->edge_vec_[edge_index];
            } else {
                assert(false);
            }
        }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        // HW0: YOUR CODE HERE
        return edge_vec_.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        if (i < num_edges()) {
            size_type e_id = edge_vec_[i].e_id;
            return Edge(this, e_id);
        } else {
            assert(false);
        }
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE
        size_type a_id = node_vec_[a.index()].n_id;
        size_type b_id = node_vec_[b.index()].n_id;
        if (a.index() >= num_nodes() || b.index()>=num_nodes()) {
            std::cout << "Accessing an incorrect node" << std::endl;
            assert(false);
        }
        
        std::pair<size_type, size_type> pair1(a_id, b_id);
        std::pair<size_type, size_type> pair2(b_id, a_id);
        bool check1 = pairs_.count(pair1);
        bool check2 = pairs_.count(pair2);
        return (check1 || check2);
    }

    /** Add an edge to the graph, or return the current edge if it already exists.
     * @pre @a a and @a b are distinct valid nodes of this graph
     * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
     * @post has_edge(@a a, @a b) == true
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
     *       Else,                        new num_edges() == old num_edges() + 1.
     *
     * Can invalidate edge indexes -- in other words, old edge(@a i) might not
     * equal new edge(@a i). Must not invalidate outstanding Edge objects.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge add_edge(const Node& a, const Node& b,
                  const edge_value_type& val = edge_value_type()) {
        // HW0: YOUR CODE HERE
        size_type a_id = node_vec_[a.index()].n_id;
        size_type b_id = node_vec_[b.index()].n_id;
        if (has_edge(a,b)) {
            //if we already have the edge, determine if we have it in orientation
            // (a,b) or orientation (b,a)
            std::pair<size_type, size_type> pair1(a_id, b_id);
            std::pair<size_type, size_type> pair2(b_id, a_id);

            // check pairs_ to see if pair1 or pair2 are keys
            bool check1 = pairs_.count(pair1);
            bool check2 = pairs_.count(pair2);
            if (check1) {
                size_type id = pairs_.at(pair1); // this is a unique edge id
                size_type idx = edge_index_map_.at(id);
                
                edge(idx).update_edge(a_id, b_id);
                return edge(idx);
            } else if (check2) {
                size_type id = pairs_.at(pair2);
                size_type idx = edge_index_map_.at(id);
                
                pairs_.erase(pair2);
                pairs_.insert({pair1, id});
                edge(idx).update_edge(a_id, b_id);
                return edge(idx);
            } else {
                assert(false);
            }
        }
        // if we don't have the edge already, determine its unique id and idx 
        size_type edge_id = edge_counter;
        edge_counter++;
        size_type edge_index = edge_vec_.size();

        // add the id -> idx mapping
        edge_index_map_.insert({edge_id, edge_index});
        
        // populate the internal edge data 
        internal_edge new_internal;
        new_internal.an_edge_ = {a_id, b_id};
        new_internal.e_id = edge_id;
        new_internal.idx_ = edge_index;
        new_internal.e_val = val;
        edge_vec_.push_back(new_internal);
        
        // add edge info to pairs_ 
        pairs_.insert({{a_id, b_id}, edge_id});

        // add b to a's neighbors and vice versa using unique id's 
        node_vec_[a.index()].neighbors.push_back(b_id);
        node_vec_[b.index()].neighbors.push_back(a_id);
        return Edge(this, edge_id);
    }
    
    /** Remove an edge from the graph, if it exists
     * @pre @a a and @a b are distinct valid nodes of this graph
     * @return 0 if the edge doesn't exist, 1 if we remove edge from graph
     * @post has_edge(@a a, @a b) == false
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
     *       Else,                        new num_edges() == old num_edges().
     *
     *
     * Complexity: O(num_nodes)
     */
    size_type remove_edge(const Node& a, const Node& b) {
        assert(has_node(a));
        assert(has_node(b));
        if (!has_edge(a,b)){ 
            //std::cout << std::endl << "there actually isn't an edge between " << node_vec_[a.index()].n_id << " and " << node_vec_[b.index()].n_id << std::endl;
            return 0;
        } else {
            size_type edge_id;
            size_type edge_idx;

            // set the unique id and idx for nodes a and b
            size_type a_idx = a.index();
            size_type a_id = node_vec_[a_idx].n_id;
            size_type b_idx = b.index();
            size_type b_id = node_vec_[b_idx].n_id;

            std::pair<size_type, size_type> pair1(a_id, b_id);        
            std::pair<size_type, size_type> pair2(b_id, a_id);     
            // check pairs_ to see if pair1 or pair2 are keys
            bool check1 = pairs_.count(pair1);

            // determine the edge id from pairs_
            if (check1) {
                edge_id = pairs_.at(pair1);
                edge_idx = edge_index_map_.at(edge_id);
                edge(edge_idx).update_edge(a_id, b_id);
            } else {
                edge_id = pairs_.at(pair2);
                edge_idx = edge_index_map_.at(edge_id);
                // if we're calling a backwards edge, switch the orientation 
                //    for the rest of the computation
                pairs_.erase(pair2);
                pairs_.insert({pair1, edge_id});
                edge(edge_idx).update_edge(a_id,b_id);
            }

            // determine the neighbors of a and the neighbors of b
            std::vector<size_type>* a_neighbor_ptr = &node_vec_[a_idx].neighbors;
            std::vector<size_type>* b_neighbor_ptr = &node_vec_[b_idx].neighbors;

            // find where b is in a and vice versa
            size_type b_remove_ = (*a_neighbor_ptr).size();
            size_type a_remove_ = (*b_neighbor_ptr).size();

            for(size_type ii=0; ii<(*a_neighbor_ptr).size(); ++ii){
                if((*a_neighbor_ptr)[ii] == b_id) {
                    b_remove_ = ii;
                }
            }
            assert(b_remove_ != (*a_neighbor_ptr).size());
            
            for(size_type ii=0; ii<(*b_neighbor_ptr).size(); ++ii){
                if((*b_neighbor_ptr)[ii] == a_id) {
                    a_remove_ = ii;
                }
            }
            assert(a_remove_ != (*b_neighbor_ptr).size());

            // swap and pop to remove b from a's neighbors
            node_vec_[a_idx].neighbors[b_remove_] = node_vec_[a_idx].neighbors.back();
            node_vec_[a_idx].neighbors.pop_back();
            
            // swap and pop to remove a from b's neighbors
            node_vec_[b_idx].neighbors[a_remove_] = node_vec_[b_idx].neighbors.back();
            node_vec_[b_idx].neighbors.pop_back();


            // overwrite edge_vec to remove edge_idx's info
            edge_vec_[edge_idx] = edge_vec_.back();
            edge_vec_[edge_idx].idx_ = edge_idx;
            // update mapping 
            edge_index_map_.at(edge_vec_[edge_idx].e_id) = edge_idx;
            
            // remove mapping for deleted edge and remove redundant edge_vec
            edge_index_map_.erase(edge_id);
            edge_vec_.pop_back();
            // remove from pairs_ because we no longer have this edge
            pairs_.erase(pair1);
            return 1;
        }
    }

    /** Remove an edge from the graph, if it exists
     * @pre @a a and @a b are distinct valid nodes of this graph
     * @return 0 if the edge doesn't exist, 1 if we remove edge from graph
     * @post has_edge(@a a, @a b) == false
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
     *       Else,                        new num_edges() == old num_edges().
     *
     *
     * Complexity: O(num_nodes)
     */
    edge_iterator remove_edge(edge_iterator e_it) {
        // dereference iterator and call other version of remove_edge 
        Edge e = *e_it;
        Node a = e.node1();
        Node b = e.node2();
        size_type x = remove_edge(a,b);
        if (x == 0) {
            return e_it;
        } else {
            return edge_begin();
        }
        
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        node_counter = 0;
        edge_counter = 0;
        node_vec_.clear();
        edge_vec_.clear();
        pairs_.clear();
        node_index_map_.clear();
        edge_index_map_.clear();
    }

    //
    // Node Iterator
    //

    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator : private totally_ordered<NodeIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node*;                    // Pointers to elements
        using reference         = Node&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::forward_iterator_tag;
                //std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid NodeIterator. */
        NodeIterator() {
        }

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:

        /** Dereference node iterator: return a Node object */
        Node operator*() const{
            if (!at_end && itr_index_<(*ptr_).size()) {
                internal_node node_info = (*ptr_)[itr_index_];
                return parent_graph_->node(node_info.idx_); 
            } else {
                return Node();
            }
        }

        /** Increment node iterator. */
        NodeIterator& operator++(){
            if (!at_end && itr_index_ < (*ptr_).size()-1) {
                itr_index_++;
                return *this;
            } else {
                itr_index_ = (*ptr_).size();
                at_end = true;
                return *this;
            }
        }

        /** Test equality of two node iterators */
        bool operator==(const NodeIterator& a) const {
            bool c1 = a.ptr_ == ptr_;
            bool c2 = a.itr_index_ == itr_index_;
            return (c1 && c2);
        }

    private:
        friend class Graph;

        // HW1 #2: YOUR CODE HERE
        graph_type* parent_graph_;
        size_type itr_index_;
        std::vector<internal_node>* ptr_;
        bool at_end = false;

        /** Constructor for initializing node iterator. */
        NodeIterator(const graph_type* parent_graph,
                     size_type itr_index,
                     const std::vector<internal_node>* ptr
                     ) :
                     parent_graph_(const_cast<graph_type*>(parent_graph)),
                     itr_index_(itr_index),
                     ptr_(const_cast<std::vector<internal_node>*>(ptr)) {
            if (itr_index == (*ptr).size()){
                at_end = true;
            }
        }
    };

    // Supply definitions AND SPECIFICATIONS for:
    /** Return node iterator starting at the beginning of node vector */
    node_iterator node_begin() const {
        return NodeIterator(this, 0, &node_vec_);
    }

    /** Return node iterator at end of node vector */
    node_iterator node_end() const {
        return NodeIterator(this, node_vec_.size(), &node_vec_);
    }

    //
    // Incident Iterator
    //

    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator : private totally_ordered<IncidentIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator() {
        }

        // HW1 #3: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        /** Dereference incident iterator returning edge between two nodes. */
        Edge operator*() const {
            if (!at_end && itr_index_ < (*ptr_).size()) {
                size_type neighbor_id = (*ptr_)[itr_index_];
                size_type curr_id = parent_graph_->node_vec_[curr_node_idx_].n_id;

                std::pair<size_type,size_type> pair1 ={curr_id, neighbor_id};
                std::pair<size_type,size_type> pair2 ={neighbor_id, curr_id};
                size_type e_id;
                if(parent_graph_->pairs_.count(pair1) > 0) {
                    e_id = parent_graph_->pairs_.at(pair1);
                } else {
                    e_id = parent_graph_->pairs_.at(pair2);
                    parent_graph_->pairs_.erase(pair2);
                    parent_graph_->pairs_.insert({pair1, e_id});
                }
                size_type e_idx = parent_graph_->edge_index_map_.at(e_id);
                parent_graph_->edge(e_idx).update_edge(curr_id, neighbor_id);
                return parent_graph_->edge(e_idx);
            } else {
                return Edge();
            }
         }

         /** Increment incident iterator to next node in neighbor list. */
        IncidentIterator& operator++(){
            if (!at_end && itr_index_ < (*ptr_).size()-1) {
                itr_index_++;
                return *this;
            } else {
                itr_index_ = (*ptr_).size();
                at_end = true;
                return *this;
            }
        }

        /** Test equality of incident iterators. */
        bool operator==(const IncidentIterator& a) const {
            bool c1 = a.ptr_ == ptr_;
            bool c2 = a.itr_index_ == itr_index_;
            bool c3 = a.curr_node_idx_ == curr_node_idx_;
            return (c1 && c2 && c3);
        }

    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE

        graph_type* parent_graph_;
        size_type curr_node_idx_;
        size_type itr_index_;
        std::vector<size_type>* ptr_;
        bool at_end = false;

        /** Private constructor for IncidentIterator. */
        IncidentIterator(const graph_type* parent_graph,
                         size_type curr_node_idx,
                         size_type itr_index,
                         const std::vector<size_type>* ptr) :
                     parent_graph_(const_cast<graph_type*>(parent_graph)),
                     curr_node_idx_(curr_node_idx),
                     itr_index_(itr_index),
                     ptr_(const_cast<std::vector<size_type>*>(ptr)) {
            if (itr_index == (*ptr).size()) {
                at_end = true;
            }
        }
    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator : private totally_ordered<EdgeIterator>{
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {
        }

        // HW1 #5: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        /** Dereference for edge iterator. Returns current Edge object. */
        Edge operator*() const {
            if (!at_end && itr_index_ < (*ptr_).size()) {
                internal_edge e = (*ptr_)[itr_index_];
                return parent_graph_->edge(e.idx_);
                //return (*ptr_).at(itr_index_);
            } else {
                return Edge();
            }
        }

        /** Increment operator for edge iterator. */
        EdgeIterator& operator++(){
            if (!at_end && itr_index_ < (*ptr_).size()-1) {
                itr_index_++;
                return *this;
            } else {
                itr_index_ = (*ptr_).size();
                at_end = true;
                return *this;
            }
        }

        /** Testing equality of edge iterators. */
        bool operator==(const EdgeIterator& a) const {
            bool c1 = a.ptr_ == ptr_;
            bool c2 = a.itr_index_ == itr_index_;
            return (c1 && c2);
        }

    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE

        graph_type* parent_graph_;
        size_type itr_index_;
        std::vector<internal_edge>* ptr_;
        bool at_end = false;

        /** Private constructor for edge iterator. */
        EdgeIterator(const graph_type* parent_graph,
                     size_type itr_index,
                     const std::vector<internal_edge>* ptr) :
                 parent_graph_(const_cast<graph_type*>(parent_graph)),
                 itr_index_(itr_index),
                 ptr_(const_cast<std::vector<internal_edge>*>(ptr)) {
            if (itr_index == (*ptr).size()) {
                at_end = true;
            }
        }
     };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return edge iterator starting at beginning of edges. */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0, &edge_vec_);
    }

    /** Return edge iterator at end of edges. */
    edge_iterator edge_end() const {
        return EdgeIterator(this, edge_vec_.size(), &edge_vec_);
    }

private:

    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.

    // Define internal node structure
    struct internal_node {
        Point pos;      // Position of a node
        size_type n_id; // Unique id for node
        size_type idx_; // Node index (between 0 and num_nodes())
        node_value_type v_type;
        std::vector<size_type> neighbors;
    };

    // Define internal edge structur e
    struct internal_edge {
        std::pair<size_type, size_type> an_edge_; //pair of unique ids
        edge_value_type e_val;
        size_type e_id;
        size_type idx_;
    };

    // vector of nodes
    std::vector<internal_node> node_vec_;
    size_type node_counter;

    // map of edges
    std::vector<internal_edge> edge_vec_;
    size_type edge_counter;

    // map of index to id
    std::map<size_type, size_type> edge_index_map_;
    std::map<size_type, size_type> node_index_map_;

    // map of pairs
    std::map<std::pair<size_type,size_type>, size_type> pairs_;

    

};

#endif // CME212_GRAPH_HPP