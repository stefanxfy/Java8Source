# Java8Source

[TOC]

## 介绍

java8源码学习注释

## 笔记整理

JUC（java.util.concurrent）是Java源码中非常重要的一个版块，无论是CAS乐观锁还是Lock悲观锁，线程池、并发集合、阻塞队列等，在日常开发中都经常用到。

如果只是停留在简单使用层面，不去深究其原理，出现了BUG，也会茫然，手足无措；而阅读其源码，了解并研究其实现原理，JUC也就不会再像一个黑盒子，平时使用也会得心应手，同时还能学习到作者的编程思维。

学习JUC源码有这么多好处，百利而无一害，何不就此开始呢？我愿与你一起学习和探讨，排除万难，领略作者匠心思维。

### atomic源码系列

疯狂撰写中...

### AQS源码系列
[AQS源码解读（番外篇）——四种自旋锁原理详解（Java代码实现SpinLock、TicketSpinLock、CLH、MCS）](https://stefan.blog.csdn.net/article/details/108750554)

[AQS源码解读（一）——AQS是什么？CLH变种体现在哪里？并发控制的核心在哪里？](https://stefan.blog.csdn.net/article/details/108817678)

[AQS源码解读（二）——从acquireQueued探索独占锁实现原理，如何阻塞？如何唤醒？](https://stefan.blog.csdn.net/article/details/108859583)

[AQS源码解读（三）——ReentrantLock原理详解（Sync、NonfairSync、FairSync）](https://stefan.blog.csdn.net/article/details/108934089)

[AQS源码解读（四）——Condition原理详解（Object#wait/notify优化？singnal唤醒线程了吗？）](https://stefan.blog.csdn.net/article/details/108946045)

[AQS源码解读（五）——从acquireShared探索共享锁实现原理，何为共享？如何共享？](https://stefan.blog.csdn.net/article/details/108950086)

[AQS源码解读（六）——从PROPAGATE和setHeadAndPropagate()分析共享锁的传播性](https://stefan.blog.csdn.net/article/details/108642253)

[AQS源码解读（七）——ReentrantReadWriteLock原理详解（读写锁是一把锁吗？如何一把锁两个状态？）](https://stefan.blog.csdn.net/article/details/108955795)

[AQS源码解读（八）——CountDownLatch倒数器原理详解](https://stefan.blog.csdn.net/article/details/108957975)

[AQS源码解读（九）——Semaphore信号量原理详解](https://stefan.blog.csdn.net/article/details/108958006)

### 并发集合

疯狂撰写中...


### 阻塞队列

疯狂撰写中...


### 线程池

疯狂撰写中...

