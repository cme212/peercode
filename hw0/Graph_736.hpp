#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
private:
    
    // HW0: YOUR CODE HERE
    // Use this space for declarations of important internal types you need
    // later in the Graph's definition.
    // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
    // code here. Just use the space if you need it.)
    
    // Declare the internal structs:
    struct internal_node;
    struct internal_edge;
    struct linked_e_n; //linked node and edge to the node
    
    
public:
    
    //
    // PUBLIC TYPE DEFINITIONS
    //
    
    /** Type of this graph. */
    using graph_type = Graph;
    
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
        // Initialize number of nodes and edges in the graph to 0.
        num_nodes_ = 0;
        num_edges_ = 0;
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
    class Node {
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
        }
        
        /** Return this node's position. */
        const Point& position() const {
            // HW0: YOUR CODE HERE
            //pointer -> vec<internal_node> --> internal_node --> Point
            return pointer_graph_-> all_nodes_[node_idx_].node_point_;
        }
        
        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            return node_idx_;
        }
        
        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // node_value_type& value();
        // const node_value_type& value() const;
        // size_type degree() const;
        // incident_iterator edge_begin() const;
        // incident_iterator edge_end() const;
        
        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            // HW0: YOUR CODE HERE
            (void) n;          // Quiet compiler warning
            // check whether it's in same graph and has valid node's index
            return ((this->pointer_graph_ == n.pointer_graph_) && (this->node_idx_ == n.node_idx_));
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
            (void) n;           // Quiet compiler warning
            // first check graph
            if(this->pointer_graph_ < n.pointer_graph_){
                return true;
            }
            // then check if in the same graph and compare node's index
            else{
                return ((this->pointer_graph_ == n.pointer_graph_) && (this->node_idx_ < n.node_idx_));
            }
        }
        
    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects
        
        // Node's pointer back to the Graph
        graph_type* pointer_graph_;
        // Node's index
        size_type node_idx_;
        // Private constructor of Node
        Node(const graph_type* g, size_type nidx){
            pointer_graph_ = const_cast<graph_type*>(g);
            node_idx_ = nidx;
        }
    };
    
    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        // HW0: YOUR CODE HERE
        return num_nodes_;
    }
    
    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }
    
    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position) {
        // HW0: YOUR CODE HERE
        
        /* step1: initialize a vector<linked_e_n> to store all nodes that linked to this new node later
         step2: create an internal_node and add it to the vector that store all nodes
         step3: update the number of nodes in the graph*/
        (void) position;      // Quiet compiler warning
        size_type idx = num_nodes_;
        std::vector<linked_e_n> linked{};
        internal_node new_node = internal_node(position, idx, linked);
        all_nodes_.push_back(new_node);
        num_nodes_ += 1;
        return Node(this, idx);
    }
    
    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        (void) n;            // Quiet compiler warning
        //check both graph and node's index
        return (this == n.pointer_graph_ && n.node_idx_ < num_nodes_);
    }
    
    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
        (void) i;             // Quiet compiler warning
        assert(i < num_nodes());
        assert(i>=0);
        return Node(this,i);
        
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
    class Edge {
    public:
        /** Construct an invalid Edge. */
        Edge() {
            // HW0: YOUR CODE HERE
        }
        
        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            if (pointer_graph_ == nullptr){
                return Node(); // Invalid Node
            }
            return Node(pointer_graph_, node1_idx_);
        }
        
        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            if (pointer_graph_ == nullptr){
                return Node(); // Invalid Node
            }
            return Node(pointer_graph_, node2_idx_);
        }
        
        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            (void) e;           // Quiet compiler warning
            //HW0: YOUR CODE HERE
            // check if the edges have same nodes
            bool check_nodes = (node1() == e.node1() && node2() == e.node2());
            bool check_nodes_rev = (node1() == e.node2() && node2() == e.node1());
            // check if the edges are in the same graph
            bool check_graph = (pointer_graph_ == e.pointer_graph_);
            // check if the edges have same index
            bool check_edge = (edge_idx_ == e.edge_idx_);
            return (check_graph && (check_nodes || check_nodes_rev) && check_edge);
        }
        
        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            (void) e;           // Quiet compiler warning
            //HW0: YOUR CODE HERE
            if (pointer_graph_ < e.pointer_graph_){
                return true;
            }else{
                return ((pointer_graph_ == e.pointer_graph_) && (this->edge_idx_ < e.edge_idx_));
            }
        }
        
    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects
        graph_type* pointer_graph_;
        size_type edge_idx_;
        size_type node1_idx_;
        size_type node2_idx_;
        // Private constructor of Edge
        Edge(const graph_type* g, const size_type eidx, const size_type n1, const size_type n2){
            pointer_graph_ = const_cast<graph_type*>(g);
            edge_idx_ = eidx;
            node1_idx_ = n1;
            node2_idx_ = n2;
        }
    };
    
    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        // HW0: YOUR CODE HERE
        return num_edges_;
    }
    
    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        (void) i;             // Quiet compiler warning
        assert(i>=0);
        assert(i< num_edges_);
        size_type nid1 = all_edges_[i].node1_idx_;
        size_type nid2 = all_edges_[i].node2_idx_;
        return Edge(this,i, nid1, nid2);
    }
    
    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE
        (void) a; (void) b;   // Quiet compiler warning
        /* step1: find the vector that contains all linked nodes and edges to Node a
         step2: loop over the vector, check if it contains node b's index
         */
        std::vector<linked_e_n> all_linked = all_nodes_[a.node_idx_].all_linked_e_n_;
        for (size_type i=0; i< all_linked.size(); i++){
            if (all_linked[i].node_idx_ == b.node_idx_){
                return true;
            }
        }
        return false;
        
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
        (void) a, (void) b;   // Quiet compiler warning
        assert(has_node(a));
        assert(has_node(b));
        assert(!(a==b));
        
        //check if the graph contains edge(a,b) already
        if(!has_edge(a,b)){
            //add the edge to the vector all_edges_
            size_type eidx = num_edges_;
            internal_edge new_edge = internal_edge(a.node_idx_, b.node_idx_,eidx);
            all_edges_.push_back(new_edge);
            
            //update linked struct linked to node a and node b
            all_nodes_[a.node_idx_].all_linked_e_n_.push_back(linked_e_n(b.node_idx_,eidx));
            all_nodes_[b.node_idx_].all_linked_e_n_.push_back(linked_e_n(a.node_idx_,eidx));
            
            //update number of edges
            num_edges_ += 1;
            return Edge(this, eidx, a.node_idx_, b.node_idx_);
        }
        return Edge();
    }
    
    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        num_nodes_ = 0;
        num_edges_ = 0;
        all_edges_.clear();
        all_nodes_.clear();
    }
    
    //
    // Node Iterator
    //
    
    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node*;                    // Pointers to elements
        using reference         = Node&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
        
        /** Construct an invalid NodeIterator. */
        NodeIterator() {
        }
        
        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Node operator*() const
        // NodeIterator& operator++()
        // bool operator==(const NodeIterator&) const
        
    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
    };
    
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_iterator node_begin() const
    // node_iterator node_end() const
    
    //
    // Incident Iterator
    //
    
    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator {
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
        // Edge operator*() const
        // IncidentIterator& operator++()
        // bool operator==(const IncidentIterator&) const
        
    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE
    };
    
    //
    // Edge Iterator
    //
    
    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator {
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
        // Edge operator*() const
        // EdgeIterator& operator++()
        // bool operator==(const EdgeIterator&) const
        
    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
    };
    
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // edge_iterator edge_begin() const
    // edge_iterator edge_end() const
    
private:
    
    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.
    
    /*Private members for class Graph*/
    size_type num_nodes_; //number of nodes in the graph
    size_type num_edges_; //number of edges in the graph
    
    /* */
    struct linked_e_n{
        size_type node_idx_; //node that linked to
        size_type edge_idx_; //edge that used to link
        
        //constructor for linked_e_n
        linked_e_n(const size_type nidx, const size_type eidx){
            node_idx_ = nidx;
            edge_idx_ = eidx;
        }
    };
    
    struct internal_node{
        Point node_point_; //node's point
        size_type node_idx_; //node's index
        std::vector<linked_e_n> all_linked_e_n_;
        
        //constructor for internal_node
        internal_node(Point p, size_type idx, std::vector<linked_e_n> linked){
            node_point_ = p;
            node_idx_ = idx;
            all_linked_e_n_ = linked;
        }
    };
    
    struct internal_edge{
        size_type node1_idx_; //node1's index
        size_type node2_idx_; //node2's index
        size_type edge_idx_; //edge's index
        
        //constructor for internal_edge
        internal_edge(const size_type n1, const size_type n2, const size_type eidx){
            node1_idx_ = n1;
            node2_idx_ = n2;
            edge_idx_ = eidx;
        }
    };
    
    // use vector<internal_nodes> & vector<internal_edges> to store all nodes and edges
    std::vector<internal_node> all_nodes_;
    std::vector<internal_edge> all_edges_;
    
};

#endif // CME212_GRAPH_HPP

