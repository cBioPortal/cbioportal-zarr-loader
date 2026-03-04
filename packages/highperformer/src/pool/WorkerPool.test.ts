import { describe, it, expect, beforeEach, vi } from 'vitest'
import { WorkerPool } from './WorkerPool'

class MockWorker {
  onmessage: ((e: MessageEvent) => void) | null = null
  postMessage = vi.fn()
  terminate = vi.fn()

  // Test helper: simulate the worker responding
  respond(data: unknown) {
    this.onmessage?.({ data } as MessageEvent)
  }
}

describe('WorkerPool', () => {
  let workers: MockWorker[]
  let pool: WorkerPool

  beforeEach(() => {
    workers = []
    const factory = () => {
      const w = new MockWorker()
      workers.push(w)
      return w as unknown as Worker
    }
    pool = new WorkerPool(factory, 2)
  })

  it('creates the requested number of workers', () => {
    expect(workers).toHaveLength(2)
  })

  it('dispatch sends message to an idle worker and resolves on response', async () => {
    const promise = pool.dispatch({ type: 'buildDefault', taskId: 'ignored' })

    // Pool assigns its own taskId — grab it from the postMessage call
    expect(workers[0].postMessage).toHaveBeenCalledTimes(1)
    const sent = workers[0].postMessage.mock.calls[0][0]
    expect(sent.type).toBe('buildDefault')
    expect(typeof sent._poolTaskId).toBe('number')

    // Simulate worker response
    workers[0].respond({ _poolTaskId: sent._poolTaskId, type: 'colorBuffer', buffer: new Uint8Array(4), version: 1 })

    const result = await promise
    expect(result.type).toBe('colorBuffer')
  })

  it('dispatches to different idle workers', async () => {
    const p1 = pool.dispatch({ type: 'buildDefault' })
    const p2 = pool.dispatch({ type: 'buildDefault' })

    expect(workers[0].postMessage).toHaveBeenCalledTimes(1)
    expect(workers[1].postMessage).toHaveBeenCalledTimes(1)

    const sent0 = workers[0].postMessage.mock.calls[0][0]
    const sent1 = workers[1].postMessage.mock.calls[0][0]

    workers[0].respond({ _poolTaskId: sent0._poolTaskId, type: 'colorBuffer', buffer: new Uint8Array(4), version: 1 })
    workers[1].respond({ _poolTaskId: sent1._poolTaskId, type: 'colorBuffer', buffer: new Uint8Array(4), version: 1 })

    await Promise.all([p1, p2])
  })

  it('queues tasks when all workers are busy and processes them when a worker frees up', async () => {
    // Fill both workers
    const p1 = pool.dispatch({ type: 'task1' })
    const p2 = pool.dispatch({ type: 'task2' })
    // This one should queue
    const p3 = pool.dispatch({ type: 'task3' })

    expect(workers[0].postMessage).toHaveBeenCalledTimes(1)
    expect(workers[1].postMessage).toHaveBeenCalledTimes(1)

    // Complete worker 0's task — should trigger queued task3
    const sent0 = workers[0].postMessage.mock.calls[0][0]
    workers[0].respond({ _poolTaskId: sent0._poolTaskId, type: 'result1' })
    await p1

    // Worker 0 should now have received the queued task
    expect(workers[0].postMessage).toHaveBeenCalledTimes(2)
    const sent0b = workers[0].postMessage.mock.calls[1][0]

    // Complete remaining tasks
    const sent1 = workers[1].postMessage.mock.calls[0][0]
    workers[1].respond({ _poolTaskId: sent1._poolTaskId, type: 'result2' })
    workers[0].respond({ _poolTaskId: sent0b._poolTaskId, type: 'result3' })

    await Promise.all([p2, p3])
  })

  it('dispose terminates all workers', () => {
    pool.dispose()
    for (const w of workers) {
      expect(w.terminate).toHaveBeenCalled()
    }
  })

  it('strips _poolTaskId from the resolved result', async () => {
    const promise = pool.dispatch({ type: 'buildDefault' })
    const sent = workers[0].postMessage.mock.calls[0][0]

    workers[0].respond({ _poolTaskId: sent._poolTaskId, type: 'colorBuffer', buffer: new Uint8Array(4), version: 1 })

    const result = await promise
    expect(result).not.toHaveProperty('_poolTaskId')
    expect(result.type).toBe('colorBuffer')
  })
})
