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
#include <string>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = int> // default template parameter
class Graph {

    public:
        /** Node value type. */
        using node_value_type = V;
        /** Edge value type. */
        using edge_value_type = E;
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
        Graph(): num_edges_ { }, num_nodes_ { }, total_nodes_ { }, nodes_map_ { },
                 adj_list_ { }, adj_map_ { }, edge_val_ { } { }

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
                Node(): uid { }, graph { } {}

                /** Return this node's position. */
                const Point& position() const {
                    // Get the position of the point with uid
                    return graph->nodes_map_[index()].position;
                }

                /** Return this node's position. */
                Point& position() {
                    // Get the position of the point with uid
                    return graph->nodes_map_[index()].position;
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    return graph->uid_idx.at(uid);
                }

                /** Return the number of nodes this node has edges connecting to. */
                size_type degree() const {
                    // Vector of indices node connects to
                    std::vector<size_type> connects = graph->adj_map_.at(uid);
                    // Size of vector
                    return (size_type)connects.size();
                }

                /** Return an iterator pointing to the first node this node
                 *  is connected to. */
                incident_iterator edge_begin() const {
                    return IncidentIterator(graph, uid, (size_type)0);
                }

                /** Return an iterator pointing to one past the last node this node
                 *  is connected to. */
                incident_iterator edge_end() const {
                    return IncidentIterator(graph, uid, degree());
                }

                /** Return this node's value */
                node_value_type& value() {
                    // Reference to node's internal structure
                    internal_node& n = const_cast<graph_type*>(graph)->nodes_map_[index()];
                    node_value_type& g = n.value;
                    return g;
                }

                /** Constant reference to node's value */
                const node_value_type& value() const {
                    return graph->nodes_map_[index()].value;
                }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    // Checks that they have the same graph and index
                    bool same_graph = n.graph == graph;
                    bool same_uid = n.uid == uid;
                    return same_graph && same_uid;
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
                    if (uid < n.uid) {
                        return true;
                    }
                    else {
                        return false;
                    }
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;
                // Uid of node
                size_type uid;
                // Pointer to graph node is apart of
                graph_type* graph;
                // Private constructor
                Node(size_type uid, const graph_type* gr): uid {uid}, graph{const_cast<graph_type*>(gr)} { }

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
            // Generate uid for new node
            size_type uid = total_nodes_;
            // This new nodes index = num_nodes_
            internal_node curr_node = internal_node {uid, position, value};
            // Map new node at index = num_nodes_ to internal structure
            nodes_map_.push_back(curr_node);
            // Map uid to empty vector in adj_map_
            adj_map_[uid] = std::vector<size_type> {};
            // Map uid to index
            uid_idx[uid] = num_nodes_;
            // Increment num_nodes_, total_nodes_
            num_nodes_++;
            total_nodes_++;
            return Node(total_nodes_ - 1, this);
        }

        /** Determine if a Node belongs to this Graph
        * @return True if @a n is currently a Node of this Graph
        *
        * Complexity: O(1).
        */
        bool has_node(const node_type& n) const {
            // Check that the graph pointed to in Node n is the current graph
            bool same_graph = this == n.graph;
            // Also check that Node n's uid is in adj_map_
            bool in_graph = adj_map_.count(n.uid) != 0;
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
            internal_node n_info = nodes_map_[i];
            node_type n = Node(n_info.uid, this);
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
                Edge(): node1_uid { }, node2_uid { }, graph { } {
                }
                /** Return a node of this Edge */
                node_type node1() const {
                    return Node(node1_uid, graph);
                }
                /** Return the other node of this Edge */
                node_type node2() const {
                    return Node(node2_uid, graph);
                }
                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const edge_type& e) const {
                    // Check associated graphs are the same
                    bool same_graph = e.graph == graph;
                    bool same_nodes = (e.node1_uid == node1_uid && e.node2_uid == node2_uid) \
                                   || (e.node1_uid == node2_uid && e.node2_uid == node1_uid);
                    return same_graph && same_nodes;
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const edge_type& e) const {
                    // Base inequality on the memory addresses of the edges
                    // Return true if memory address is smaller than memory
                    // address of e
                    unsigned long curr_add = (unsigned long)&(*this);
                    unsigned long e_add = (unsigned long)&e;
                    return curr_add < e_add;
                }

                /** Return length of edge */
                double length() const {
                    return norm(node1().position() - node2().position());
                }

                /** Return value associated with edge */
                edge_value_type& value() {
                    // Smaller uid maps to larger uid
                    size_type smaller = std::min(node1_uid, node2_uid);
                    size_type larger = std::max(node1_uid, node2_uid);
                    edge_value_type& edge_info = graph->edge_val_[smaller][larger];
                    return edge_info;
                }

                /** Return value associated with edge */
                const edge_value_type& value() const {
                    // Smaller uid maps to larger uid
                    size_type smaller = std::min(node1_uid, node2_uid);
                    size_type larger = std::max(node1_uid, node2_uid);
                    edge_value_type& edge_info = graph->edge_val_[smaller][larger];
                    return edge_info;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;
                // Edge goes from node uid node1_uid to node2_uid
                size_type node1_uid;
                size_type node2_uid;
                // Graph associated with edge
                graph_type* graph;
                // Private constructor
                Edge(const graph_type* gr, size_type uid1, size_type uid2):
                     node1_uid {uid1}, node2_uid {uid2}, graph {const_cast<graph_type*>(gr)} { }

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
            int n1_uid = adj_list_[i][0];
            int n2_uid = adj_list_[i][1];
            return Edge(this, n1_uid, n2_uid);
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

            // Return false if uid has no neighboring nodes
            if (adj_map_.at(a.uid).size() == 0) {
                return false;
            }

            // Loop through Node a's neighbors and return true if Node b's
            // index is found
            std::vector<size_type> neighbors = adj_map_.at(a.uid);
            for (size_type idx = 0; idx < neighbors.size(); idx++) {
                if (neighbors[idx] == b.uid) {
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
                // Add tuple containing connected uids to adj_list_
                adj_list_.emplace_back(std::vector<size_type>{a.uid, b.uid});
                // Add to adj_map_
                adj_map_[a.uid].emplace_back(b.uid);
                adj_map_[b.uid].emplace_back(a.uid);
                // Add to edge_val_
                edge_val_[std::min(a.uid, b.uid)][std::max(a.uid, b.uid)] = edge_value_type{};
                num_edges_ += 1;
            }
            return Edge(this, a.uid, b.uid);
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
            uid_idx.clear();
            edge_val_.clear();
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
                graph_type* graph;
                // Index of node
                size_type node_idx;
                // Private constructor
                NodeIterator(const graph_type* gr, size_type idx) :
                    graph {const_cast<graph_type*>(gr)}, node_idx {idx} {}
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
                    return Edge(graph, root_uid, connects[con_idx]);
                }

                /** IncidentIterators are equal if they have the same associated
                 *  graphs, the same root node uid, and the current edge
                 *  leads to the same node */
                bool operator==(const IncidentIterator& edge) const {
                    return (edge.graph == graph && edge.root_uid == root_uid &&
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
                graph_type* graph;
                // uid of root node in associated graph
                size_type root_uid;
                // Vector of uids of nodes the root node is connected to
                std::vector<size_type> connects;
                // Index in neighboring nodes
                size_type con_idx;
                // Number of nodes root node is connected to
                size_type num_connects;

                // Private constructor
                IncidentIterator(const graph_type* gr, size_type uid, size_type con_idx) :
                graph{const_cast<graph_type*>(gr)}, root_uid{uid}, con_idx{con_idx} {
                    // num_connects is degree of root node
                    num_connects = gr->node(gr->uid_idx.at(root_uid)).degree();
                    // If root node has neighbors, get vector of their indices
                    if (num_connects != 0) {
                        connects = gr->adj_map_.at(root_uid);
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
                    /** Iterators equal if they have the same associated graph
                      * and edge index */
                    return (ei.graph == graph && ei.edge_idx == edge_idx);
                }

            private:
                friend class Graph;
                // Graph associated with edges
                graph_type* graph;
                // Index of current edge
                size_type edge_idx;
                // Private constructor
                EdgeIterator(const graph_type* gr, size_type idx) :
                    graph{const_cast<graph_type*>(gr)}, edge_idx{idx} {}
        };

        /** Return EdgeIterator pointing to the first edge in the graph */
        edge_iterator edge_begin() const {
            return EdgeIterator(this, 0);
        }

        /** Return EdgeIterator pointing to one past the last edge in the graph */
        edge_iterator edge_end() const {
            return EdgeIterator(this, num_edges());
        }

        /**
         * @brief Removes the edge between two nodes
         *
         * @param[in] a One node the edge is connected to
         * @param[in] b The other node the edge is connected to
         * @return Status code of function run
         *
         * @post has_edge(a, b) == False
         * @post If old has_edge(a, b) == False, then old num_edges() == num_edges()
         * @post If old has_edge(a, b) == True, then old num_edges() == num_edges() + 1
         *
         * Complexity: O(num_edges + degree(a) + degree(b))
         *
         */
        size_type remove_edge(const Node& a, const Node& b) {
            // Return if there isn't an edge between nodes
            if (!has_edge(a, b)) {
                return 0;
            }
            // Find and remove edge in adj_list_
            for (unsigned int i = 0; i < adj_list_.size(); i++) {
                if ((adj_list_[i][0] == a.uid && adj_list_[i][1] == b.uid) ||
                    (adj_list_[i][1] == a.uid && adj_list_[i][0] == b.uid)) {
                        adj_list_.erase(adj_list_.begin() + i);
                        break;
                    }
            }
            // Remove from adj_map_
            for (unsigned int i = 0; i < adj_map_[a.uid].size(); i++) {
                if (adj_map_[a.uid][i] == b.uid) {
                    adj_map_[a.uid].erase(adj_map_[a.uid].begin() + i);
                    break;
                }
            }
            // Remove "symmetric" entries from adj_map_ as well
            for (unsigned int i = 0; i < adj_map_[b.uid].size(); i++) {
                if (adj_map_[b.uid][i] == a.uid) {
                    adj_map_[b.uid].erase(adj_map_[b.uid].begin() + i);
                    break;
                }
            }
            // Decrement number of edges
            num_edges_--;
            return 1;
        }

        /**
         * @brief Removes the edge between two nodes
         *
         * @param[in] e The edge to remove
         * @return Status code of function run
         *
         * @post has_edge(e.node1(), e.node2()) == False
         * @post If old has_edge(e.node1(), e.node2()) == False, then old num_edges() == num_edges()
         * @post If old has_edge(e.node1(), e.node2()) == True, then old num_edges() == num_edges() + 1
         *
         * Complexity: O(num_edges + degree(e.node1()) + degree(e.node2()))
         */
        size_type remove_edge(const Edge& e) {
            return remove_edge(e.node1(), e.node2());
        }

        /**
         * @brief Removes the edge between two nodes
         *
         * @param[in] e_it An edge_iterator pointing at the edge _e_ to remove
         * @return Status code of function run
         *
         * @post has_edge((_e_).node1(), (_e_).node2()) == False
         * @post If old has_edge((_e_).node1(), (_e_).node2()) == False, then old num_edges() == num_edges()
         * @post If old has_edge((_e_).node1(), (_e_).node2()) == True, then old num_edges() == num_edges() + 1
         *
         * Complexity: O(num_edges + degree(_e_.node1()) + degree(_e_.node2()))
         */
        edge_iterator remove_edge(edge_iterator e_it) {
            // Dereference and remove edge iterator is pointing to
            remove_edge(*e_it);
            return e_it;
        }

        /**
         * @brief Removes a node in the graph
         *
         * @param[in] n Node to remove
         * @return Status code of function run
         *
         * @post has_node(_n_) == False
         * @post If old has_node(_n_) == False, then old num_nodes() == num_nodes()
         * @post If old has_node(_n_) == True, then old num_edges() == num_edges() + 1
         *
         * Complexity: O(degree(_n_) * (num_edges + degree(_n_) + degree(neighbor)))
         * Removing only nodes while keeping edges has constant time complexity.
         * There are degree(_n_) edges, and removing each edge has time complexity
         * O(num_edges + degree(a) + degree(neighbor)), where there is an edge from
         * node a to a neighboring node.
         */
        size_type remove_node(const Node& n) {
            // Return if node isn't in graph
            if (!has_node(n)) {
                return 0;
            }

            // Get current index of node before it is removed
            size_type node_idx = n.index();

            // Remove edge between n and incident nodes
            std::vector<size_type> connects = adj_map_.at(n.uid);

            for (unsigned int i = 0; i < connects.size(); i++) {
                node_type Node2 = Node(connects[i], this);
                remove_edge(n, Node2);
            }

            // Remove node from nodes_map_
            nodes_map_.erase(nodes_map_.begin() + uid_idx.at(n.uid));
            // Decrement number of nodes
            num_nodes_--;

            // Update other nodes' indices in uid_idx map
            for (size_type idx = node_idx; idx < num_nodes_; idx++) {
                size_type curr_uid = nodes_map_[idx].uid;
                uid_idx[curr_uid] = idx;
            }

            return 1;
        }

        /**
         * @brief Removes a node in the graph
         *
         * @param[in] n_it node_iterator pointing to the node _n_ to remove
         * @return Status code of function run
         *
         * @post has_node(_n_) == False
         * @post If old has_node(_n_) == False, then old num_nodes() == num_nodes()
         * @post If old has_node(_n_) == True, then old num_edges() == num_edges() + 1
         *
         * Complexity: O(degree(_n_) * (num_edges + degree(_n_) + degree(neighbor)))
         */
        node_iterator remove_node(node_iterator n_it ) {
            // Dereference and remove node node iterator is pointing to
            remove_node(*n_it);
            return n_it;
        }

    private:

        // Number of edges and nodes in a graph
        size_type num_edges_;
        size_type num_nodes_;
        // Total nodes. The total number of nodes added over the graph's lifetime
        // Used to generate node uids (as re-indexing doens't effect node uid,
        // removing nodes doesn't effect total_nodes_)
        size_type total_nodes_;

        // Internal details of a node.
        struct internal_node {
            size_type uid;
            Point position;
            node_value_type value;
        };

        // Maps node uid to corresponding internal_node
        // Index in the vector is the index in the graph
        std::vector<internal_node> nodes_map_;
        // Adjacency list. adj_list_[i] will have a vector consisting of the
        // uids of the two nodes connected by edge i
        std::vector<std::vector<size_type>> adj_list_;
        // Map node uid to uids of other nodes it connects to
        std::unordered_map<size_type, std::vector<size_type>> adj_map_;
        // Map node uid to the node's index in the graph
        std::unordered_map<size_type, size_type> uid_idx;
        // Map a node uid to another map, mapping connected node uids to the edge value
        std::unordered_map<size_type, std::unordered_map<size_type, edge_value_type>> edge_val_;
};

#endif // CME212_GRAPH_HPP
