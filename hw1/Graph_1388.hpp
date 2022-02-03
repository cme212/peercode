#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int> // default template parameter
class Graph {

    public:
        /** Node value type. */
        using node_value_type = V;
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
        Graph(): num_edges_ { 0 }, num_nodes_ { 0 }, nodes_map_ { },
                 adj_list_ { }, adj_map_ { } { }

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
        class Node : private totally_ordered<Node> {
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
                Node(): idx { }, graph { } {}

                /** Return this node's position. */
                const Point& position() const {
                    // Get the position of the point with idx
                    return graph->nodes_map_.at(idx).position;
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    return idx;
                }

                /** Return the number of nodes this node has edges connecting to. */
                size_type degree() const {
                    // If node's index is not in adj_map_ then it has no neighbors */
                    if (graph->adj_map_.count(index()) == 0) {
                        return (size_type)0;
                    }
                    // Vector of indices node connects to
                    std::vector<size_type> connects = graph->adj_map_.at(index());
                    // Size of vector
                    return (size_type)connects.size();
                }

                /** Return an iterator pointing to the first node this node
                 *  is connected to. */
                incident_iterator edge_begin() const {
                    return IncidentIterator(graph, idx, (size_type)0);
                }

                /** Return an iterator pointing to one past the last node this node
                 *  is connected to. */
                incident_iterator edge_end() const {
                    return IncidentIterator(graph, idx, degree());
                }

                /** Return this node's value */
                node_value_type& value() {
                    // Reference to node's internal structure
                    internal_node& n = const_cast<graph_type*>(graph)->nodes_map_.at(idx);
                    node_value_type& g = n.value;
                    return g;
                }

                /** Constant reference to node's value */
                const node_value_type& value() const {
                    return graph->nodes_map_.at(idx).value;
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    // Checks that they have the same graph and index
                    bool same_graph = n.graph == graph;
                    bool same_index = n.index() == index();
                    return same_graph && same_index;
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
                    // Compare their indexes
                    if (idx < n.index()) {
                        return true;
                    }
                    // Same index but different graph, consider node with the
                    // smaller graph
                    if (idx == n.index() && n.graph != graph) {
                        return graph->num_nodes() < n.graph->num_nodes();
                    }
                    else {
                        return false;
                    }
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;
                // Index of node in graph
                size_type idx;
                // Pointer to graph node is apart of
                const graph_type* graph;
                // Private constructor
                Node(size_type idx, const graph_type* gr): idx {idx}, graph {gr} { }

            };

        /** Return the number of nodes in the graph.
        *
        * Complexity: O(1).
        */
        size_type size() const {
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
        node_type add_node(const Point& position, const node_value_type& value = node_value_type()) {
            // This new nodes index = num_nodes_
            internal_node curr_node = internal_node {num_nodes_, position, value};
            // Map new node at index = num_nodes_ to internal structure
            nodes_map_[num_nodes_] = curr_node;
            // Increment num_nodes_
            num_nodes_+=1;
            return Node(num_nodes_ - 1, this);
        }

        /** Determine if a Node belongs to this Graph
        * @return True if @a n is currently a Node of this Graph
        *
        * Complexity: O(1).
        */
        bool has_node(const node_type& n) const {
            // Check that the graph pointed to in Node n is the current graph
            bool same_graph = this == n.graph;
            // Also check that Node n's index is in the current graph
            bool in_graph = nodes_map_.count(n.idx) != 0;
            return same_graph && in_graph;
        }

        /** Return the node with index @a i.
        * @pre 0 <= @a i < num_nodes()
        * @post result_node.index() == i
        *
        * Complexity: O(1).
        */
        node_type node(size_type i) const {
            assert(i < num_nodes());
            // Get information of node at index i
            internal_node n_info = nodes_map_.at(i);
            node_type n =  Node(n_info.index, this);

            assert(n.index() == i);

            return n;
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
        class Edge : private totally_ordered<Edge> {
            public:
                /** Construct an invalid Edge. */
                Edge(): node1_idx { }, node2_idx { }, graph { } {
                }
                /** Return a node of this Edge */
                node_type node1() const {
                    return Node(node1_idx, graph);
                }
                /** Return the other node of this Edge */
                node_type node2() const {
                    return Node(node2_idx, graph);
                }
                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const edge_type& e) const {
                    // Check associated graphs are the same
                    bool same_graph = e.graph == graph;
                    bool same_nodes = (e.node1_idx == node1_idx && e.node2_idx == node2_idx) \
                                   || (e.node1_idx == node2_idx && e.node2_idx == node1_idx);
                    return same_graph && same_nodes;
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const edge_type& e) const {
                    // Base inequality on the sum of the indices of the nodes making
                    // up the edges
                    size_type idx_sum1 = node1_idx + node2_idx;
                    size_type idx_sum2 = e.node1_idx + e.node2_idx;
                    if (idx_sum1 < idx_sum2) {
                        return true;
                    }
                    return false;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;
                // Edge goes from node index node1_idx to node2_idx
                size_type node1_idx;
                size_type node2_idx;
                // Graph associated with edge
                const graph_type* graph;
                // Private constructor
                Edge(const graph_type* gr, size_type idx1, size_type idx2):
                     node1_idx {idx1}, node2_idx {idx2}, graph {gr} { }

            };

        /** Return the total number of edges in the graph.
        *
        * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
        */
        size_type num_edges() const {
            return num_edges_;
        }

        /** Return the edge with index @a i.
        * @pre 0 <= @a i < num_edges()
        *
        * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
        */
        edge_type edge(size_type i) const {
            assert(i < num_edges());
            int n1_idx = std::get<0>(adj_list_[i]);
            int n2_idx = std::get<1>(adj_list_[i]);
            return Edge(this, n1_idx, n2_idx);
        }

        /** Test whether two nodes are connected by an edge.
        * @pre @a a and @a b are valid nodes of this graph
        * @return True if for some @a i, edge(@a i) connects @a a and @a b.
        *
        * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
        */
        bool has_edge(const node_type& a, const node_type& b) const {

            // Return false if the graph doesn't contain both nodes
            if (!has_node(a) || !has_node(b)) {
                return false;
            }

            // Return false of index has no neighboring nodes
            if (adj_map_.count(a.index()) == 0) {
                return false;
            }

            // Loop through Node a's neighbors and return true if Node b's
            // index is found
            std::vector<size_type> neighbors = adj_map_.at(a.index());
            for (size_type idx = 0; idx < neighbors.size(); idx++) {
                if (neighbors[idx] == b.index()) {
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
        edge_type add_edge(const node_type& a, const node_type& b) {
            // Return invalid edge if both nodes aren't in the graph
            if (!has_node(a) || !has_node(b)) {
                return Edge();
            }
            // Add edge if it doesn't exist
            if (!has_edge(a, b)) {
                // Add tuple containing connected indices to adj_list_
                adj_list_.emplace_back(std::tuple<size_type, size_type>{a.idx, b.idx});
                // If there is no edge between node a and any other node
                if (adj_map_.count(a.idx) == 0) {
                    adj_map_[a.idx] = std::vector<size_type>();
                }
                adj_map_[a.idx].emplace_back(b.idx);
                // If there is no edge between node b and any other node
                if (adj_map_.count(b.idx) == 0) {
                    adj_map_[b.idx] = std::vector<size_type>();
                }
                adj_map_[b.idx].emplace_back(a.idx);
                num_edges_+=1;
            }
            return Edge(this, a.idx, b.idx);
        }

        /** Remove all nodes and edges from this graph.
        * @post num_nodes() == 0 && num_edges() == 0
        *
        * Invalidates all outstanding Node and Edge objects.
        */
        void clear() {
            // Clear or zero out defining attributes of the graph
            adj_list_.clear();
            nodes_map_.clear();
            adj_map_.clear();
            num_edges_ = 0;
            num_nodes_ = 0;
        }

        //
        // Node Iterator
        //

        /** @class Graph::NodeIterator
        * @brief Iterator class for nodes. A forward iterator. */
        class NodeIterator : private totally_ordered<NodeIterator>{
            public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Node;                     // Element type
                using pointer           = Node*;                    // Pointers to elements
                using reference         = Node&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy
                /** Construct an invalid NodeIterator. */
                NodeIterator() : graph{}, node_idx{} {}

                /** Dereference iterator to return Node */
                value_type operator*() const {
                    return graph->node(node_idx);
                }

                /** Define ++ operator for iterator */
                NodeIterator& operator++() {
                    // increment index
                    node_idx++;
                    return *this;
                }

                /** NodeIterators equal if they point to the same graph and
                 * index */
                bool operator==(const NodeIterator& iter) const {
                    bool graph_same = iter.graph == graph;
                    bool idx_same = iter.node_idx == node_idx;
                    return graph_same && idx_same;
                }

            private:
                friend class Graph;
                // Associated graph of node
                const graph_type* graph;
                // Index of node
                size_type node_idx;
                // Private constructor
                NodeIterator(const graph_type* gr, size_type idx) : graph{gr}, node_idx{idx} {}
        };

        /** Return NodeIterator pointing to the first node in the graph */
        node_iterator node_begin() const {
            return NodeIterator(this, 0);
        }

        /** Return NodeIterator pointing to one past the last node in the graph */
        node_iterator node_end() const {
            return NodeIterator(this, num_nodes_);
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
                using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid IncidentIterator. */
                IncidentIterator() {}

                /** Dereference iterator and return the associated edge */
                value_type operator*() const {
                    // Construct edge such that root node is the starting node
                    return Edge(graph, root_idx, connects[con_idx]);
                }

                /** IncidentIterators are equal if they have the same associated
                 *  graphs, the same root node index, and the current edge
                 *  leads to the same node */
                bool operator==(const IncidentIterator& edge) const {
                    return (edge.graph == graph && edge.root_idx == root_idx &&
                            edge.con_idx == con_idx);
                }

                /** Define ++ operator for IncidentIterator */
                IncidentIterator& operator++() {
                    // Increment connected index
                    con_idx++;
                    return *this;
                }

            private:
                friend class Graph;
                // Associated graph
                const graph_type* graph;
                // Index of root node in associated graph
                size_type root_idx;
                // Vector of indices of nodes the root node is connected to
                std::vector<size_type> connects;
                // Index in neighboring nodes
                size_type con_idx;
                // Number of nodes root node is connected to
                size_type num_connects;

                // Private constructor
                IncidentIterator(const graph_type* gr, size_type idx, size_type con_idx) :
                graph{gr}, root_idx{idx}, con_idx{con_idx} {
                    // num_connects is degree of root node
                    num_connects = gr->node(root_idx).degree();
                    // If root node has neighbors, get vector of their indices
                    if (num_connects != 0) {
                        connects = gr->adj_map_.at(root_idx);
                    }
                    else {
                        connects = {};
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
                using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid EdgeIterator. */
                EdgeIterator() : graph{ }, edge_idx{ } {}

                /** Dereference iterate and return the current edge */
                value_type operator*() const {
                    return graph->edge(edge_idx);
                }

                /** Define ++ operator to move to the next edge */
                EdgeIterator& operator++() {
                    edge_idx++;
                    return *this;
                }

                /** Define equality for EdgeIterators */
                bool operator==(const EdgeIterator& ei) const {
                    /** Iterators equal if they have the same associated grave
                      * and edge index */
                    return (ei.graph == graph && ei.edge_idx == edge_idx);
                }

            private:
                friend class Graph;
                // Graph associated with edges
                const graph_type* graph;
                // Index of current edge
                size_type edge_idx;
                // Private constructor
                EdgeIterator(const graph_type* gr, size_type idx) :
                    graph{gr}, edge_idx{idx} {}
        };

        /** Return EdgeIterator pointing to the first edge in the graph */
        edge_iterator edge_begin() const {
            return EdgeIterator(this, 0);
        }

        /** Return EdgeIterator pointing to one past the last edge in the graph */
        edge_iterator edge_end() const {
            return EdgeIterator(this, num_edges());
        }

    private:

        // Number of edges and nodes in a graph
        size_type num_edges_;
        size_type num_nodes_;

        // Internal details of a node.
        struct internal_node {
            size_type index;
            Point position;
            node_value_type value;
        };

        // Maps node index to corresponding internal_node
        std::unordered_map<size_type, internal_node> nodes_map_;
        // Adjacency list. adj_list_[i] will have a tuple consisting of the
        // indexes of the two nodes connected by edge i
        std::vector<std::tuple<size_type, size_type>> adj_list_;
        // Map node index to indexes of other nodes it connects to
        std::unordered_map<size_type, std::vector<size_type>> adj_map_;


};

#endif // CME212_GRAPH_HPP
