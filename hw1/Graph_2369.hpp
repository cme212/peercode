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
template <typename V = int>
class Graph {
private:

    // HW0: YOUR CODE HERE
    // Use this space for declarations of important internal types you need
    // later in the Graph's definition.
    // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
    // code here. Just use the space if you need it.)

    // Declare struct for node's internal data
    struct internal_node;

public:

    //
    // PUBLIC TYPE DEFINITIONS
    //

    /** Type of node value. */
    using node_value_type = V;

    /** Type of this graph. */
    using graph_type = Graph<V>;

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

        /** Return this node's position. */
        const Point& position() const {
            // HW0: YOUR CODE HERE
            return fetch().pos;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            return fetch().n_id;
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
            return IncidentIterator(parent_graph_,n_id_,0,&fetch().neighbors);
        }

        /** Return the end of node's incident iterator over neighbors. */
        incident_iterator edge_end() const {
            return IncidentIterator(parent_graph_, n_id_,
                                    fetch().neighbors.size(),
                                    &fetch().neighbors);
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            // HW0: YOUR CODE HERE
            bool c1 = n_id_ == n.index();
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
            if (this->n_id_ < n.index()) {
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
            // check if n_id_ is a valid id
            if (n_id_ < parent_graph_->node_vec_.size()) {
                return parent_graph_->node_vec_[n_id_];
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
        size_type new_id = node_counter;
        node_counter++;
        internal_node add_internal;
        add_internal.pos = position;
        add_internal.n_id = new_id;
        add_internal.v_type = val;
        add_internal.neighbors = {};
        node_vec_.push_back(add_internal);
        Node add = Node(this, new_id);
        return add;
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
            return Node(this,i);
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

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            size_type node1_idx = an_edge_.first;
            return (*parent_graph_).node(node1_idx);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            size_type node2_idx = an_edge_.second;
            return (*parent_graph_).node(node2_idx);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            // HW0: YOUR CODE HERE
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
            } else if ((node1() < e.node1()) && (node2() < e.node2())) {
                return true;
            }
            return false;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members sand methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects

        // An edge has a pair of node indices
        std::pair<size_type, size_type> an_edge_;
        // An edge has a pointer to the graph it is in
        Graph* parent_graph_;
        // An edge has an id associated with it
        size_type e_id_;
        // Flag to see if Edge is valid
        bool is_valid_;

        /** Private constructor (used by Graph) creates valid edge. */
        Edge(graph_type* parent_graph, const Node& a,
                                       const Node& b, size_type e_id) {
            parent_graph_ = parent_graph;
            an_edge_.first = a.index();
            an_edge_.second = b.index();
            e_id_ = e_id;
            is_valid_ = true;
        }

        /** Helper method to update the current orientation of edge. */
        void update_edge(size_type a, size_type b) {
            an_edge_.first = a;
            an_edge_.second = b;
        }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        // HW0: YOUR CODE HERE
        return edges_.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        return edges_.at(i);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE
        std::pair<size_type, size_type> pair1(a.index(), b.index());
        std::pair<size_type, size_type> pair2(b.index(), a.index());
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
    Edge add_edge(const Node& a, const Node& b) {
        // HW0: YOUR CODE HERE
        if (has_edge(a,b)) {
            std::pair<size_type, size_type> pair1(a.index(), b.index());
            std::pair<size_type, size_type> pair2(b.index(), a.index());

            // check pairs_ to see if pair1 or pair2 are keys
            bool check1 = pairs_.count(pair1);
            bool check2 = pairs_.count(pair2);
            if (check1) {
                size_type id = pairs_.at(pair1);
                edges_.at(id).update_edge(a.index(), b.index());
                return edges_.at(id);
            } else if (check2) {
                size_type id = pairs_.at(pair2);
                edges_.at(id).update_edge(a.index(), b.index());
                return edges_.at(id);
            } else {
                std::cout << "Something went wrong with add_edge" << std::endl;
            }
        }
        size_type edge_id = edge_counter;
        edge_counter++;

        Edge new_edge = Edge(this, a, b, edge_id);
        edges_.insert({edge_id, new_edge});
        pairs_.insert({{a.index(),b.index()}, edge_id});
        // add b to a's neighbors and vice versa
        node_vec_.at(a.index()).neighbors.push_back(b.index());
        node_vec_.at(b.index()).neighbors.push_back(a.index());
        return new_edge;
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
        edges_.clear();
        pairs_.clear();
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
                return parent_graph_->node(node_info.n_id);
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
                const Node n1 = parent_graph_->node(curr_node_id_);
                const Node n2 = parent_graph_->node(neighbor_id);
                return parent_graph_->add_edge(n1, n2);
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
            bool c3 = a.curr_node_id_ == curr_node_id_;
            return (c1 && c2 && c3);
        }

    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE

        graph_type* parent_graph_;
        size_type curr_node_id_;
        size_type itr_index_;
        std::vector<size_type>* ptr_;
        bool at_end = false;

        /** Private constructor for IncidentIterator. */
        IncidentIterator(const graph_type* parent_graph,
                         size_type curr_node_id,
                         size_type itr_index,
                         const std::vector<size_type>* ptr) :
                     parent_graph_(const_cast<graph_type*>(parent_graph)),
                     curr_node_id_(curr_node_id),
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
                return (*ptr_).at(itr_index_);
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
        std::map<size_type, Edge>* ptr_;
        bool at_end = false;

        /** Private constructor for edge iterator. */
        EdgeIterator(const graph_type* parent_graph,
                     size_type itr_index,
                     const std::map<size_type, Edge>* ptr) :
                 parent_graph_(const_cast<graph_type*>(parent_graph)),
                 itr_index_(itr_index),
                 ptr_(const_cast<std::map<size_type,Edge>*>(ptr)) {
            if (itr_index == (*ptr).size()) {
                at_end = true;
            }
        }
     };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return edge iterator starting at beginning of edges. */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0, &edges_);
    }

    /** Return edge iterator at end of edges. */
    edge_iterator edge_end() const {
        return EdgeIterator(this, edges_.size(), &edges_);
    }

private:

    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.

    // Define internal node structure
    struct internal_node {
        Point pos;      // Position of a node
        size_type n_id; // Unique id for node
        node_value_type v_type;
        std::vector<size_type> neighbors;
    };

    // map of nodes
    std::vector<internal_node> node_vec_;
    size_type node_counter;

    // map of edges
    std::map<size_type, Edge> edges_;
    size_type edge_counter;

    // map of pairs
    std::map<std::pair<size_type,size_type>, size_type> pairs_;
};

#endif // CME212_GRAPH_HPP
